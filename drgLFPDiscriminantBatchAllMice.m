function drgLFPDiscriminantBatchAllMice

%Discriminant and PCA analysis for LFP data
%Performs one LDA analysis using data from all mice

% It takes as input a choices file such as
% drgbChoicesDiscriminantJustin_spm_perfom_LFP_20180215


% handles.drgbchoices.which_discriminant chooses the analysis:
%
%1 Perceptron for power LFP (very slow and has not been troublehsot)
%
%2 Linear discriminant analysis (LDA) for power LFP
%
%3 Principal component analysis (PCA) for power LFP
%
%4 Linear discriminant analysis (LDA) for phase in phase amplitude
%       coupling (PAC)
%
%5 Principal component analysis for PAC
%
%6 LDA for subsets of electrodes for power LFP
%
%7 LDA for power calculated with the high freq Hilbert at the PAC phase
%
%8 PCA for power calculated with the high freq Hilbert at the PAC phase
%
%9 LDA for power calculated with the high freq Hilbert at the PAC phase
%
%10 LDA for wavelet power calculated at the PAC phase
%
%11 PCA for wavelet power calculated at the PAC phase
%
%12 LDA for subsets of electrodes for wavelet power calculated at the PAC phase
%
%13 LDA for no laser/laser for wavelet power calculated at the PAC phase
%
%14 PCA for no laser/laser for wavelet power calculated at the PAC phase
%
%15 LDA for wavelet power calculated at the PAC phase using data for all
%mice
%
%16 PCA for wavelet power calculated at the PAC phase using data for all
%mice
%
%17 LDA for subsets of electrodes for wavelet power calculated at the PAC
%phase for all mice (not vetted yet)

close all
clear all

percent_correct_ii=1;

first_file=1;
figNo=0;

[choiceFileName,choiceBatchPathName] = uigetfile({'drgbChoicesDiscriminant*.m'},'Select the .m file with all the choices for discriminant analysis');
fprintf(1, ['\ndrgLFPDiscriminantBatch run for ' choiceFileName '\n\n']);

tempDirName=['temp' choiceFileName(12:end-2)];

addpath(choiceBatchPathName)
eval(['handles=' choiceFileName(1:end-2) ';'])
handles.choiceFileName=choiceFileName;
handles.choiceBatchPathName=choiceBatchPathName;

p=gcp;

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

%Uncomment to find out which files were processed
for mouseNo=1:max(handles.drgbchoices.mouse_no)
    for groupNo=1:max(handles.drgbchoices.group_no)
        
        %Find out whether this mouse/group was processed
        if (sum(handles.drgbchoices.which_discriminant==10)>0)
            try
                dcalc=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo);
                process_data_for_mouse=0;
                fprintf(1, ['Mouse no %d ' handles.drgbchoices.group_no_names{groupNo} ' was processed\n\n'],mouseNo);
            catch
                process_data_for_mouse=1;
                fprintf(1, ['Mouse no %d ' handles.drgbchoices.group_no_names{groupNo} ' was not processed\n\n'],mouseNo);
            end
        end
    end
end

%Did the user specify keep_phase?
if ~isfield(handles.drgbchoices,'keep_phase')
    handles.drgbchoices.keep_phase=0;
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

phase_per_file=[];

if all_files_present==1
    
    all_mouse_wav=[];
    
    handles_out.drgbchoices=handles.drgbchoices;
    
    handles.displayData=0;
    
    handles.subtractRef=handles.drgbchoices.subtractRef;
    handles.evTypeNo=handles.drgbchoices.referenceEvent; %Process all OdorOn trials
    
    no_files=handles.drgbchoices.no_files;
    
    for mouseNo=1:max(handles.drgbchoices.mouse_no)
        if (sum(handles.drgbchoices.which_discriminant==13)>0)
            no_trialsPACwave=0;
            mouse_has_data=0;
        end
        
        for groupNo=1:max(handles.drgbchoices.group_no)
            
            
            fprintf(1, ['Mouse no %d ' handles.drgbchoices.group_no_names{groupNo} '\n\n'],mouseNo);
            
            process_data_for_mouse=1;
            if overwrite_out_file==1
                process_data_for_mouse=1;
            else
                %Find out whether this mouse/group was processed
                if (sum(handles.drgbchoices.which_discriminant==2)>0)
                    try
                        dcalc=handles_out.discriminant_per_mouse(mouseNo).group(groupNo);
                        process_data_for_mouse=0;
                    catch
                        process_data_for_mouse=1;
                    end
                end
                %Find out whether this mouse/group was processed
                if (sum(handles.drgbchoices.which_discriminant==10)>0)
                    try
                        dcalc=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo);
                        process_data_for_mouse=0;
                    catch
                        process_data_for_mouse=1;
                    end
                end
                %Find out whether this mouse/group was processed
                if (sum(handles.drgbchoices.which_discriminant==13)>0)
                    try
                        dcalc=handles_out.discriminant_PACpower_per_mouse(mouseNo).group(1).which_analysis(1).discriminant_calculated;
                        process_data_for_mouse=0;
                    catch
                        process_data_for_mouse=1;
                    end
                end
            end
            
            if process_data_for_mouse==1
                no_trials=0;
                no_trialsPAC=0;
                
                if ~(sum(handles.drgbchoices.which_discriminant==13)>0)
                    no_trialsPACwave=0;
                    mouse_has_data=0;
                end
                
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
                            jtPathNameSource=handles.drgbchoices.PathName{filNum};
                        else
                            jtPathNameSource=handles.drgbchoices.PathName;
                        end
                        
                        if exist([jtPathNameSource jtFileName])==0
                            waitfor(msgbox('Please plug the hard drive back in...'));
                        end
                        
                        jtPathName=handles.drgb.outPathName;
                        
                        %Copy the file to the ouput directory
                        copyfile([jtPathNameSource jtFileName],[jtPathName jtFileName]);
                        
                        this_rhd=[jtFileName(10:end-4) '.rhd'];
                        if exist([jtPathNameSource this_rhd])==0
                            this_dg=[jtFileName(10:end-4) '.dg'];
                            copyfile([jtPathNameSource this_dg],[jtPathName this_dg]);
                            while exist([jtPathName this_dg])==0
                            end
                            %                                 handles.drg.drta_p.fullName=[jtPathName jtFileName(10:end-4) '.dg'];
                        else
                            this_rhd=[jtFileName(10:end-4) '.rhd'];
                            copyfile([jtPathNameSource this_rhd],[jtPathName this_rhd]);
                            while exist([jtPathName this_rhd])==0
                            end
                            %                                 handles.drg.drta_p.fullName=[jtPathName jtFileName(10:end-4) '.rhd'];
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
                        %Initalize the output variable for parfor for LFP
                        %power
                        if (sum(handles.drgbchoices.which_discriminant==1)>0)||(sum(handles.drgbchoices.which_discriminant==2)>0)...
                                ||(sum(handles.drgbchoices.which_discriminant==3)>0)||(sum(handles.drgbchoices.which_discriminant==6)>0)
                            
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
                        end
                        
                        %Initalize the output variable for parfor for PAC
                        if (sum(handles.drgbchoices.which_discriminant==4)>0)||(sum(handles.drgbchoices.which_discriminant==5)>0)...
                                ||(sum(handles.drgbchoices.which_discriminant==7)>0)||(sum(handles.drgbchoices.which_discriminant==10)>0)...
                                ||(sum(handles.drgbchoices.which_discriminant==13)>0)||(sum(handles.drgbchoices.which_discriminant==15)>0)...
                                ||(sum(handles.drgbchoices.which_discriminant==16)>0)
                            for LFPNo=1:no_elect
                                
                                for ii=1:handles.drgbchoices.no_PACpeaks
                                    %Enter the per LFP values
                                    par_out(LFPNo).PAC(ii).no_trials=0;
                                    par_out(LFPNo).PAC(ii).meanVectorLength=[];
                                    par_out(LFPNo).PAC(ii).meanVectorAngle=[];
                                    par_out(LFPNo).PAC(ii).peakAngle=[];
                                    par_out(LFPNo).PAC(ii).mod_indx=[];
                                    par_out(LFPNo).PAC(ii).this_trialNo=[];
                                    par_out(LFPNo).PAC(ii).all_phase_histo=[];
                                    par_out(LFPNo).PAC(ii).perCorrPAC=[];
                                    par_out(LFPNo).PAC(ii).which_eventPAC=[];
                                    par_out(LFPNo).PAC(ii).meanPeakAngle=[];
                                    par_out(LFPNo).PAC(ii).PACtimecourse=[];
                                    par_out(LFPNo).PAC(ii).peakAngleForPower=[];
                                    par_out(LFPNo).PAC(ii).troughAngleForPower=[];
                                    par_out(LFPNo).t_power=[];
                                end
                            end
                            
                        end
                        
                        if (sum(handles.drgbchoices.which_discriminant==10)>0)||(sum(handles.drgbchoices.which_discriminant==13)>0)...
                                ||(sum(handles.drgbchoices.which_discriminant==15)>0)
                            for LFPNo=1:no_elect
                                
                                for ii=1:handles.drgbchoices.no_PACpeaks
                                    %Enter the per LFP values
                                    par_out(LFPNo).t_pac=[];
                                    par_out(LFPNo).PACwave(ii).PACtimecourse=[];
                                    par_out(LFPNo).PACwave(ii).meanPeakPower=[];
                                    par_out(LFPNo).PACwave(ii).meanTroughPower=[];
                                    par_out(LFPNo).PACwave(ii).meanPower=[];
                                    
                                end
                            end
                            
                        end
                        
                        
                        parfor LFPNo=1:no_elect
                            %                         for LFPNo=1:no_elect
                            handlespf=struct();
                            handlespf=handles;
                            
                            %Now calculate the LFP power
                            if (sum(handlespf.drgbchoices.which_discriminant==1)>0)||(sum(handlespf.drgbchoices.which_discriminant==2)>0)...
                                    ||(sum(handlespf.drgbchoices.which_discriminant==3)>0)||(sum(handlespf.drgbchoices.which_discriminant==6)>0)
                                this_peakLFPNo=handlespf.drgbchoices.which_electrodes(LFPNo);
                                handlespf.peakLFPNo=this_peakLFPNo;
                                handlespf.burstLFPNo=this_peakLFPNo;
                                handlespf.lastTrialNo=handlespf.drg.session(1).noTrials;
                                all_Power=[];
                                all_Power_ref=[];
                                all_Power_timecourse=[];
                                perCorr_pertr=[];
                                which_event=[];
                                [t,f,all_Power,all_Power_ref, all_Power_timecourse, this_trialNo, perCorr_pertr, which_event]=drgGetLFPPowerForThisEvTypeNo(handlespf);
                                
                                %Note: This prunning of trials is here because for a small number
                                %of sessions trials included are nether S+ nor S-
                                %e.g in M4_spmc_170518_095813 only 83 of 119 are OK
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
                                
                                
                                
                                fprintf(1, 'For file no %d, LFP no %d the number of trials included for LFP Power is %d (out of %d)\n',filNum,this_peakLFPNo,length(this_trialNo),handlespf.lastTrialNo);
                                
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
                            
                            %Now calculate the PAC
                            if (sum(handlespf.drgbchoices.which_discriminant==4)>0)||(sum(handlespf.drgbchoices.which_discriminant==5)>0)||...
                                    (sum(handlespf.drgbchoices.which_discriminant==7)>0)||(sum(handlespf.drgbchoices.which_discriminant==8)>0)||...
                                    (sum(handlespf.drgbchoices.which_discriminant==9)>0)
                                for PACii=1:handlespf.drgbchoices.no_PACpeaks
                                    this_peakLFPNo=handlespf.drgbchoices.which_electrodes(LFPNo);
                                    handlespf.peakLFPNo=this_peakLFPNo;
                                    handlespf.burstLFPNo=this_peakLFPNo;
                                    handlespf.n_peak=PACii;
                                    handlespf.peakLowF=handlespf.drgbchoices.PACpeakLowF;
                                    handlespf.peakHighF=handlespf.drgbchoices.PACpeakHighF;
                                    handlespf.burstLowF=handlespf.drgbchoices.PACburstLowF(PACii);
                                    handlespf.burstHighF=handlespf.drgbchoices.PACburstHighF(PACii);
                                    handlespf.n_phase_bins=handlespf.drgbchoices.n_phase_bins;
                                    handlespf.lastTrialNo=handlespf.drg.session(1).noTrials;
                                    
                                    %Please note this is the same function called by
                                    %drgMaster when the user chooses Phase Amplitude
                                    %Coupling
                                    handlespf=drgThetaAmpPhaseTrialRange(handlespf);
                                    
                                    %Note: This prunning of trials is here because for a small number
                                    %of sessions trials included are nether S+ nor S-
                                    %e.g in M4_spmc_170518_095813 only 83 of 119 are OK
                                    this_trialNo=handlespf.drgb.PAC.this_trialNo;
                                    which_event=handlespf.drgb.PAC.which_event;
                                    these_all_which_events=zeros(length(handles.drgbchoices.events_to_discriminate),length(this_trialNo));
                                    for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                        kk=find(handles.drgbchoices.evTypeNos==handles.drgbchoices.events_to_discriminate(ii));
                                        these_all_which_events(ii,:)= which_event(kk,:);
                                    end
                                    
                                    valid_trials=logical(sum(these_all_which_events,1));
                                    
                                    old_meanVectorLength=handlespf.drgb.PAC.meanVectorLength;
                                    meanVectorLength=zeros(1,sum(valid_trials));
                                    meanVectorLength(:,:)=old_meanVectorLength(1,valid_trials);
                                    
                                    old_meanVectorAngle=handlespf.drgb.PAC.meanVectorAngle;
                                    meanVectorAngle=zeros(1,sum(valid_trials));
                                    meanVectorAngle(:,:)=old_meanVectorAngle(1,valid_trials);
                                    
                                    old_peakAngle=handlespf.drgb.PAC.peakAngle;
                                    peakAngle=zeros(1,sum(valid_trials));
                                    peakAngle(:,:)=old_peakAngle(1,valid_trials);
                                    
                                    old_mod_indx=handlespf.drgb.PAC.mod_indx;
                                    mod_indx=zeros(1,sum(valid_trials));
                                    mod_indx(:,:)=old_mod_indx(1,valid_trials);
                                    
                                    old_this_trialNo=handlespf.drgb.PAC.this_trialNo;
                                    this_trialNo=zeros(1,sum(valid_trials));
                                    this_trialNo(:,:)=old_this_trialNo(1,valid_trials);
                                    
                                    old_all_phase_histo=handlespf.drgb.PAC.all_phase_histo;
                                    szaph=size(old_all_phase_histo);
                                    all_phase_histo=zeros(sum(valid_trials),szaph(2));
                                    all_phase_histo(:,:)=old_all_phase_histo(valid_trials,:);
                                    
                                    old_perCorr=handlespf.drgb.PAC.perCorr;
                                    perCorr=zeros(1,sum(valid_trials));
                                    perCorr(:,:)=old_perCorr(1,valid_trials);
                                    
                                    old_which_event=handlespf.drgb.PAC.which_event;
                                    szwe=size(old_which_event);
                                    which_event=zeros(szwe(1),sum(valid_trials));
                                    which_event(:,:)=old_which_event(:,valid_trials);
                                    
                                    old_meanPeakAngle=handlespf.drgb.PAC.meanPeakAngle;
                                    meanPeakAngle=zeros(1,sum(valid_trials));
                                    meanPeakAngle(:,:)=old_meanPeakAngle(1,valid_trials);
                                    
                                    old_meanPeakPower=handlespf.drgb.PAC.meanPeakPower;
                                    meanPeakPower=zeros(1,sum(valid_trials));
                                    meanPeakPower(:,:)=old_meanPeakPower(1,valid_trials);
                                    
                                    old_meanTroughPower=handlespf.drgb.PAC.meanTroughPower;
                                    meanTroughPower=zeros(1,sum(valid_trials));
                                    meanTroughPower(:,:)=old_meanTroughPower(1,valid_trials);
                                    
                                    jj=0;
                                    out_PACtimecourse=struct;
                                    for kk=1:length(handlespf.drgb.PAC.PACtimecourse)
                                        if valid_trials(kk)==1
                                            jj=jj+1;
                                            out_PACtimecourse(jj).out_times=handlespf.drgb.PAC.PACtimecourse(kk).out_times;
                                            out_PACtimecourse(jj).out_phase=handlespf.drgb.PAC.PACtimecourse(kk).out_phase;
                                            out_PACtimecourse(jj).out_time_PAChisto=handlespf.drgb.PAC.PACtimecourse(kk).out_time_PAChisto;
                                            out_PACtimecourse(jj).peakPower=handlespf.drgb.PAC.PACtimecourse(kk).peakPower;
                                            out_PACtimecourse(jj).troughPower=handlespf.drgb.PAC.PACtimecourse(kk).troughPower;
                                        end
                                    end
                                    
                                    %Enter the per LFP values
                                    par_out(LFPNo).PAC(PACii).no_trials=sum(valid_trials);
                                    par_out(LFPNo).PAC(PACii).meanVectorLength=meanVectorLength;
                                    par_out(LFPNo).PAC(PACii).meanVectorAngle=meanVectorAngle;
                                    par_out(LFPNo).PAC(PACii).peakAngle=peakAngle;
                                    par_out(LFPNo).PAC(PACii).mod_indx=mod_indx;
                                    par_out(LFPNo).PAC(PACii).this_trialNo=this_trialNo;
                                    par_out(LFPNo).PAC(PACii).all_phase_histo=all_phase_histo;
                                    par_out(LFPNo).PAC(PACii).perCorrPAC=perCorr;
                                    par_out(LFPNo).PAC(PACii).which_eventPAC=which_event;
                                    par_out(LFPNo).PAC(PACii).meanPeakAngle=meanPeakAngle;
                                    par_out(LFPNo).PAC(PACii).meanPeakPower=meanPeakPower;
                                    par_out(LFPNo).PAC(PACii).meanTroughPower=meanTroughPower;
                                    par_out(LFPNo).PAC(PACii).PACtimecourse=out_PACtimecourse;
                                    par_out(LFPNo).t_power=handlespf.drgb.PAC.t_power;
                                    
                                    
                                    
                                end
                                fprintf(1, 'For file no %d, LFP no %d the number of trials included in PAC is %d (out of %d)\n',filNum,this_peakLFPNo,par_out(LFPNo).PAC(PACii).no_trials,handlespf.lastTrialNo);
                                
                            end
                            
                            %Now calculate wave power and PAC
                            if (sum(handlespf.drgbchoices.which_discriminant==10)>0)||(sum(handlespf.drgbchoices.which_discriminant==11)>0)...
                                    ||(sum(handlespf.drgbchoices.which_discriminant==12)>0)||(sum(handles.drgbchoices.which_discriminant==13)>0)...
                                    ||(sum(handles.drgbchoices.which_discriminant==14)>0)||(sum(handles.drgbchoices.which_discriminant==15)>0)...
                                    ||(sum(handles.drgbchoices.which_discriminant==16)>0)
                                for PACii=1:handlespf.drgbchoices.no_PACpeaks
                                    this_peakLFPNo=handlespf.drgbchoices.which_electrodes(LFPNo);
                                    handlespf.peakLFPNo=this_peakLFPNo;
                                    handlespf.burstLFPNo=this_peakLFPNo;
                                    handlespf.n_peak=PACii;
                                    handlespf.peakLowF=handlespf.drgbchoices.PACpeakLowF;
                                    handlespf.peakHighF=handlespf.drgbchoices.PACpeakHighF;
                                    handlespf.burstLowF=handlespf.drgbchoices.PACburstLowF(PACii);
                                    handlespf.burstHighF=handlespf.drgbchoices.PACburstHighF(PACii);
                                    handlespf.n_phase_bins=handlespf.drgbchoices.n_phase_bins;
                                    handlespf.lastTrialNo=handlespf.drg.session(1).noTrials;
                                    
                                    if handlespf.drgbchoices.keep_phase==1
                                        if filNum~=handles.drgbchoices.referenceAngleFile(filNum)
                                            handlespf.use_peakAngle=1;
                                            refFilNum=handles.drgbchoices.referenceAngleFile(filNum);
                                            handlespf.peakAngle_for_power=phase_per_file(refFilNum).LFP(LFPNo).PAC(PACii).peakAngleForPower;
                                            handlespf.troughAngle_for_power=phase_per_file(refFilNum).LFP(LFPNo).PAC(PACii).troughAngleForPower;
                                        else
                                            handlespf.use_peakAngle=0;
                                        end
                                    end
                                    
                                    %Please note this is the same function called by
                                    %drgMaster when the user chooses Phase Amplitude
                                    %Coupling
                                    handlespf=drgLFPwaveTimecourse(handlespf);
                                    
                                    if handlespf.drgb.PAC.no_trials>0
                                        
                                        %Note: This prunning of trials is here because for a small number
                                        %of sessions trials included are nether S+ nor S-
                                        %e.g in M4_spmc_170518_095813 only 83 of 119 are OK
                                        this_trialNo=handlespf.drgb.PAC.this_trialNo;
                                        which_event=handlespf.drgb.PAC.which_event;
                                        these_all_which_events=zeros(length(handles.drgbchoices.events_to_discriminate),length(this_trialNo));
                                        for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                            kk=find(handles.drgbchoices.evTypeNos==handles.drgbchoices.events_to_discriminate(ii));
                                            these_all_which_events(ii,:)= which_event(kk,:);
                                        end
                                        
                                        valid_trials=logical(sum(these_all_which_events,1));
                                        
                                        old_meanVectorLength=handlespf.drgb.PAC.meanVectorLength;
                                        meanVectorLength=zeros(1,sum(valid_trials));
                                        meanVectorLength(:,:)=old_meanVectorLength(1,valid_trials);
                                        
                                        old_meanVectorAngle=handlespf.drgb.PAC.meanVectorAngle;
                                        meanVectorAngle=zeros(1,sum(valid_trials));
                                        meanVectorAngle(:,:)=old_meanVectorAngle(1,valid_trials);
                                        
                                        old_peakAngle=handlespf.drgb.PAC.peakAngle;
                                        peakAngle=zeros(1,sum(valid_trials));
                                        peakAngle(:,:)=old_peakAngle(1,valid_trials);
                                        
                                        old_mod_indx=handlespf.drgb.PAC.mod_indx;
                                        mod_indx=zeros(1,sum(valid_trials));
                                        mod_indx(:,:)=old_mod_indx(1,valid_trials);
                                        
                                        old_this_trialNo=handlespf.drgb.PAC.this_trialNo;
                                        this_trialNo=zeros(1,sum(valid_trials));
                                        this_trialNo(:,:)=old_this_trialNo(1,valid_trials);
                                        
                                        old_all_phase_histo=handlespf.drgb.PAC.all_phase_histo;
                                        szaph=size(old_all_phase_histo);
                                        all_phase_histo=zeros(sum(valid_trials),szaph(2));
                                        all_phase_histo(:,:)=old_all_phase_histo(valid_trials,:);
                                        
                                        old_perCorr=handlespf.drgb.PAC.perCorr;
                                        perCorr=zeros(1,sum(valid_trials));
                                        perCorr(:,:)=old_perCorr(1,valid_trials);
                                        
                                        old_which_event=handlespf.drgb.PAC.which_event;
                                        szwe=size(old_which_event);
                                        which_event=zeros(szwe(1),sum(valid_trials));
                                        which_event(:,:)=old_which_event(:,valid_trials);
                                        
                                        old_meanPeakAngle=handlespf.drgb.PAC.meanPeakAngle;
                                        meanPeakAngle=zeros(1,sum(valid_trials));
                                        meanPeakAngle(:,:)=old_meanPeakAngle(1,valid_trials);
                                        
                                        old_meanPeakPower=handlespf.drgb.PAC.meanPeakPower;
                                        meanPeakPower=zeros(1,sum(valid_trials));
                                        meanPeakPower(:,:)=old_meanPeakPower(1,valid_trials);
                                        
                                        old_meanTroughPower=handlespf.drgb.PAC.meanTroughPower;
                                        meanTroughPower=zeros(1,sum(valid_trials));
                                        meanTroughPower(:,:)=old_meanTroughPower(1,valid_trials);
                                        
                                        jj=0;
                                        out_PACtimecourse=struct;
                                        for kk=1:length(handlespf.drgb.PAC.PACtimecourse)
                                            if valid_trials(kk)==1
                                                jj=jj+1;
                                                out_PACtimecourse(jj).out_times=handlespf.drgb.PAC.PACtimecourse(kk).out_times;
                                                out_PACtimecourse(jj).out_phase=handlespf.drgb.PAC.PACtimecourse(kk).out_phase;
                                                out_PACtimecourse(jj).out_time_PAChisto=handlespf.drgb.PAC.PACtimecourse(kk).out_time_PAChisto;
                                                out_PACtimecourse(jj).peakPower=handlespf.drgb.PAC.PACtimecourse(kk).peakPower;
                                                out_PACtimecourse(jj).troughPower=handlespf.drgb.PAC.PACtimecourse(kk).troughPower;
                                            end
                                        end
                                        
                                        %Now do PACwave
                                        old_meanPeakPowerwave=handlespf.drgb.PACwave.meanPeakPower;
                                        meanPeakPowerwave=zeros(1,sum(valid_trials));
                                        meanPeakPowerwave(:,:)=old_meanPeakPowerwave(1,valid_trials);
                                        
                                        old_meanTroughPowerwave=handlespf.drgb.PACwave.meanTroughPower;
                                        meanTroughPowerwave=zeros(1,sum(valid_trials));
                                        meanTroughPowerwave(:,:)=old_meanTroughPowerwave(1,valid_trials);
                                        
                                        old_meanPowerwave=handlespf.drgb.PACwave.meanPower;
                                        meanPowerwave=zeros(1,sum(valid_trials));
                                        meanPowerwave(:,:)=old_meanPowerwave(1,valid_trials);
                                        
                                        jj=0;
                                        out_PACtimecoursewave=struct;
                                        for kk=1:length(handlespf.drgb.PACwave.PACtimecourse)
                                            if valid_trials(kk)==1
                                                jj=jj+1;
                                                out_PACtimecoursewave(jj).peakPower=handlespf.drgb.PACwave.PACtimecourse(kk).peakPower;
                                                out_PACtimecoursewave(jj).troughPower=handlespf.drgb.PACwave.PACtimecourse(kk).troughPower;
                                                out_PACtimecoursewave(jj).troughPower=handlespf.drgb.PACwave.PACtimecourse(kk).troughPower;
                                            end
                                        end
                                        
                                        %Enter the per LFP values
                                        par_out(LFPNo).PAC(PACii).no_trials=sum(valid_trials);
                                        par_out(LFPNo).PAC(PACii).meanVectorLength=meanVectorLength;
                                        par_out(LFPNo).PAC(PACii).meanVectorAngle=meanVectorAngle;
                                        par_out(LFPNo).PAC(PACii).peakAngle=peakAngle;
                                        par_out(LFPNo).PAC(PACii).mod_indx=mod_indx;
                                        par_out(LFPNo).PAC(PACii).this_trialNo=this_trialNo;
                                        par_out(LFPNo).PAC(PACii).all_phase_histo=all_phase_histo;
                                        par_out(LFPNo).PAC(PACii).perCorrPAC=perCorr;
                                        par_out(LFPNo).PAC(PACii).which_eventPAC=which_event;
                                        par_out(LFPNo).PAC(PACii).meanPeakAngle=meanPeakAngle;
                                        par_out(LFPNo).PAC(PACii).meanPeakPower=meanPeakPower;
                                        par_out(LFPNo).PAC(PACii).meanTroughPower=meanTroughPower;
                                        par_out(LFPNo).PAC(PACii).PACtimecourse=out_PACtimecourse;
                                        par_out(LFPNo).t_power=handlespf.drgb.PAC.t_power;
                                        par_out(LFPNo).PAC(PACii).peakAngleForPower=handlespf.drgb.PAC.peakAngleForPower;
                                        par_out(LFPNo).PAC(PACii).troughAngleForPower=handlespf.drgb.PAC.troughAngleForPower;
                                        
                                        %Enter the per LFP values for PACwave
                                        par_out(LFPNo).PACwave(PACii).meanPeakPower=meanPeakPowerwave;
                                        par_out(LFPNo).PACwave(PACii).meanTroughPower=meanTroughPowerwave;
                                        par_out(LFPNo).PACwave(PACii).meanPower=meanPowerwave;
                                        par_out(LFPNo).PACwave(PACii).PACtimecourse=out_PACtimecoursewave;
                                        par_out(LFPNo).t_pac=handlespf.drgb.PACwave.t_pac;
                                    end
                                    
                                end
                                fprintf(1, 'For file no %d, LFP no %d the number of trials included in PAC is %d (out of %d)\n',filNum,this_peakLFPNo,par_out(LFPNo).PAC(PACii).no_trials,handlespf.lastTrialNo);
                                
                            end
                            
                        end
                        
                        t_power=par_out(1).t_power;
                        t_pac=par_out(1).t_pac;
                        
                        
                        %Extract all_log_P_timecourse for powerLFP
                        if (sum(handles.drgbchoices.which_discriminant==1)>0)||(sum(handles.drgbchoices.which_discriminant==2)>0)...
                                ||(sum(handles.drgbchoices.which_discriminant==3)>0)||(sum(handles.drgbchoices.which_discriminant==6)>0)
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
                                        szt=size(this_all_log_P_timecourse,4);
                                        all_log_P_timecourse=zeros(length(handles.drgbchoices.lowF),length(handles.drgbchoices.which_electrodes),par_out(LFPNo).length_trial_no+no_trials,length(t));
                                        all_log_P_timecourse(:,:,1:no_trials,1:szt)=this_all_log_P_timecourse(:,:,1:no_trials,:);
                                    end
                                end
                                for bwii=1:length(handles.drgbchoices.lowF)
                                    all_log_P_timecourse(bwii,LFPNo,no_trials+1:no_trials+par_out(LFPNo).length_trial_no,:)=par_out(LFPNo).bwii(bwii).log_P_timecourse;
                                end
                            end
                        end
                        
                        %Extract all_phase_timecourse for PAC phase
                        if (sum(handles.drgbchoices.which_discriminant==4)>0)||(sum(handles.drgbchoices.which_discriminant==5)>0)
                            t_pac=par_out(1).PAC(1).PACtimecourse(1).out_times;
                            no_bins=size(par_out(1).PAC(1).PACtimecourse(1).out_time_PAChisto,2);
                            for LFPNo=1:length(handles.drgbchoices.which_electrodes)
                                
                                %                             this_peakLFPNo=handles.drgbchoices.which_electrodes(LFPNo);
                                if LFPNo>1
                                    if last_No_trials~=par_out(LFPNo).PAC(1).no_trials
                                        fprintf(1, 'WARNING. For file no %d the number of trials differ between different LFPs\n ',filNum);
                                    end
                                end
                                last_No_trials=par_out(LFPNo).PAC(1).no_trials;
                                
                                if LFPNo==1
                                    if no_trialsPAC==0
                                        %all_phase_timecourse=zeros(length(handles.drgbchoices.PACburstLowF),length(handles.drgbchoices.which_electrodes),par_out(LFPNo).PAC(1).no_trials,length(t_pac));
                                        all_phase_histo_timecourse=zeros(length(handles.drgbchoices.PACburstLowF),length(handles.drgbchoices.which_electrodes)*no_bins,par_out(LFPNo).PAC(1).no_trials,length(t_pac));
                                    else
                                        %                                         this_all_phase_timecourse=[];
                                        %                                         this_all_phase_timecourse=all_phase_timecourse;
                                        %                                         all_phase_timecourse=zeros(length(handles.drgbchoices.PACburstLowF),length(handles.drgbchoices.which_electrodes),par_out(LFPNo).PAC(1).no_trials+no_trialsPAC,length(t_pac));
                                        %                                         all_phase_timecourse(:,:,1:no_trialsPAC,:)=this_all_phase_timecourse;
                                        
                                        this_all_phase_histo_timecourse=[];
                                        this_all_phase_histo_timecourse=all_phase_histo_timecourse;
                                        all_phase_histo_timecourse=zeros(length(handles.drgbchoices.PACburstLowF),length(handles.drgbchoices.which_electrodes)*no_bins,par_out(LFPNo).PAC(1).no_trials+no_trialsPAC,length(t_pac));
                                        all_phase_histo_timecourse(:,:,1:no_trialsPAC,:)=this_all_phase_histo_timecourse;
                                    end
                                end
                                for PACii=1:length(handles.drgbchoices.PACburstLowF)
                                    %                                     these_phases=zeros(par_out(LFPNo).PAC(1).no_trials,length(t_pac));
                                    %                                     for trNo=1:par_out(LFPNo).PAC(1).no_trials
                                    %                                         these_phases(trNo,:)=par_out(LFPNo).PAC(PACii).PACtimecourse(trNo).out_phase;
                                    %                                     end
                                    %                                     all_phase_timecourse(PACii,LFPNo,no_trialsPAC+1:no_trialsPAC+par_out(LFPNo).PAC(1).no_trials,:)=these_phases;
                                    
                                    for trNo=1:par_out(LFPNo).PAC(1).no_trials
                                        all_phase_histo_timecourse(PACii,(LFPNo-1)*no_bins+1:LFPNo*no_bins,no_trialsPAC+trNo,:)=par_out(LFPNo).PAC(PACii).PACtimecourse(trNo).out_time_PAChisto';
                                    end
                                    
                                end
                            end
                        end
                        
                        %Extract timecourse for PAC peak power
                        if (sum(handles.drgbchoices.which_discriminant==7)>0)||(sum(handles.drgbchoices.which_discriminant==8)>0)||...
                                (sum(handles.drgbchoices.which_discriminant==9)>0)
                            t=t_power;
                            for LFPNo=1:length(handles.drgbchoices.which_electrodes)
                                
                                %                             this_peakLFPNo=handles.drgbchoices.which_electrodes(LFPNo);
                                if LFPNo>1
                                    if last_No_trials~=length(par_out(LFPNo).PAC(1).meanPeakPower)
                                        fprintf(1, 'WARNING. For file no %d the number of trials differ between different LFPs\n ',filNum);
                                    end
                                end
                                last_No_trials=length(par_out(LFPNo).PAC(1).meanPeakPower);
                                
                                if LFPNo==1
                                    if no_trials==0
                                        all_log_P_timecoursePACpeak=zeros(length(handles.drgbchoices.PACburstLowF),length(handles.drgbchoices.which_electrodes),length(par_out(LFPNo).PAC(1).meanPeakPower),length(t));
                                        all_log_P_timecoursePACtrough=zeros(length(handles.drgbchoices.PACburstLowF),length(handles.drgbchoices.which_electrodes),length(par_out(LFPNo).PAC(1).meanTroughPower),length(t));
                                    else
                                        %Peak
                                        this_all_log_P_timecoursePAC=[];
                                        this_all_log_P_timecoursePAC=all_log_P_timecoursePACpeak;
                                        szt=size(this_all_log_P_timecoursePAC,4);
                                        all_log_P_timecoursePACpeak=zeros(length(handles.drgbchoices.PACburstLowF),length(handles.drgbchoices.which_electrodes),length(par_out(LFPNo).PAC(1).meanPeakPower)+no_trials,length(t));
                                        all_log_P_timecoursePACpeak(:,:,1:no_trials,1:szt)=this_all_log_P_timecoursePAC(:,:,1:no_trials,:);
                                        %Trough
                                        this_all_log_P_timecoursePAC=[];
                                        this_all_log_P_timecoursePAC=all_log_P_timecoursePACtrough;
                                        szt=size(this_all_log_P_timecoursePAC,4);
                                        all_log_P_timecoursePACtrough=zeros(length(handles.drgbchoices.PACburstLowF),length(handles.drgbchoices.which_electrodes),length(par_out(LFPNo).PAC(1).meanTroughPower)+no_trials,length(t));
                                        all_log_P_timecoursePACtrough(:,:,1:no_trials,1:szt)=this_all_log_P_timecoursePAC(:,:,1:no_trials,:);
                                    end
                                end
                                for PACii=1:length(handles.drgbchoices.PACburstLowF)
                                    for trNo=1:length(par_out(LFPNo).PAC(1).meanPeakPower)
                                        all_log_P_timecoursePACpeak(PACii,LFPNo,no_trials+trNo,:)=par_out(LFPNo).PAC(PACii).PACtimecourse(trNo).peakPower;
                                    end
                                    for trNo=1:length(par_out(LFPNo).PAC(1).meanTroughPower)
                                        all_log_P_timecoursePACtrough(PACii,LFPNo,no_trials+trNo,:)=par_out(LFPNo).PAC(PACii).PACtimecourse(trNo).troughPower;
                                    end
                                end
                                
                                
                            end
                        end
                        
                        %Extract timecourse for wavelet power calculated at PAC phase
                        if (sum(handles.drgbchoices.which_discriminant==10)>0)||(sum(handles.drgbchoices.which_discriminant==11)>0)||...
                                (sum(handles.drgbchoices.which_discriminant==12)>0)||(sum(handles.drgbchoices.which_discriminant==13)>0)...
                                ||(sum(handles.drgbchoices.which_discriminant==14)>0)||(sum(handles.drgbchoices.which_discriminant==15)>0)...
                                ||(sum(handles.drgbchoices.which_discriminant==16)>0)
                            
                            %Save the peak and trough phase
                            for LFPNo=1:length(handles.drgbchoices.which_electrodes)
                                for PACii=1:length(handles.drgbchoices.PACburstLowF)
                                    phase_per_file(filNum).LFP(LFPNo).PAC(PACii).peakAngleForPower=par_out(LFPNo).PAC(PACii).peakAngleForPower;
                                    phase_per_file(filNum).LFP(LFPNo).PAC(PACii).troughAngleForPower=par_out(LFPNo).PAC(PACii).troughAngleForPower;
                                end
                            end
                            
                            %If there are different numbers of trials weed
                            %out the mismatched trials
                            PACii=1;
                            lastNoTrials=par_out(1).PAC(PACii).no_trials;
                            these_trials=par_out(1).PAC(PACii).this_trialNo;
                            trials_match=1;
                            minLFP=1;
                            for LFPNo=2:length(handles.drgbchoices.which_electrodes)
                                if par_out(LFPNo).PAC(PACii).no_trials~=lastNoTrials
                                    trials_match=0;
                                end
                                these_trials=intersect(these_trials,par_out(LFPNo).PAC(PACii).this_trialNo);
                            end
                            
                            if trials_match==0
                                %The code below makes all LFPs have the
                                %same number of trials shared by lte
                                %largest subset of electrodes
                                
                                %Find largest subset of shared trials
                                PACii=1;
                                no_trn=1;
                                no_trials_per_LFP(no_trn)=par_out(1).PAC(PACii).no_trials;
                                no_LFPs(no_trn)=1;
                                for LFPNo=2:length(handles.drgbchoices.which_electrodes)
                                    if sum(no_trials_per_LFP==par_out(LFPNo).PAC(PACii).no_trials)>0
                                        ii=find(no_trials_per_LFP==par_out(LFPNo).PAC(PACii).no_trials);
                                        no_LFPs(ii)=no_LFPs(ii)+1;
                                    else
                                        no_trn=no_trn+1;
                                        no_LFPs(no_trn)=1;
                                        no_trials_per_LFP(no_trn)=par_out(LFPNo).PAC(PACii).no_trials;
                                    end
                                end
                                
                                %Find the largest subset
                                [max_no, max_ii]=max(no_LFPs);
                                maxLFPs=zeros(1,max_no);
                                jj=0;
                                for LFPNo=1:length(handles.drgbchoices.which_electrodes)
                                    if no_trials_per_LFP(max_ii)==par_out(LFPNo).PAC(PACii).no_trials
                                        jj=jj+1;
                                        maxLFPs(jj)=LFPNo;
                                    end
                                end
                                
                                for LFPNo=1:length(handles.drgbchoices.which_electrodes)
                                    for PACii=1:length(handles.drgbchoices.PACburstLowF)
                                        if  par_out(LFPNo).PAC(PACii).no_trials~=no_trials_per_LFP(max_ii)
                                            [min_dist,ii_replace]=min(abs(LFPNo-maxLFPs));
                                            LFPrep=maxLFPs(ii_replace(1));
                                            par_out(LFPNo).PAC(PACii).no_trials=par_out(LFPrep).PAC(PACii).no_trials;
                                            par_out(LFPNo).PAC(PACii).this_trialNo=par_out(LFPrep).PAC(PACii).this_trialNo;
                                            par_out(LFPNo).PAC(PACii).meanVectorLength=par_out(LFPrep).PAC(PACii).meanVectorLength;
                                            par_out(LFPNo).PAC(PACii).meanVectorAngle=par_out(LFPrep).PAC(PACii).meanVectorAngle;
                                            par_out(LFPNo).PAC(PACii).peakAngle=par_out(LFPrep).PAC(PACii).peakAngle;
                                            par_out(LFPNo).PAC(PACii).mod_indx=par_out(LFPrep).PAC(PACii).mod_indx;
                                            par_out(LFPNo).PAC(PACii).all_phase_histo=par_out(LFPrep).PAC(PACii).all_phase_histo;
                                            par_out(LFPNo).PAC(PACii).perCorrPAC=par_out(LFPrep).PAC(PACii).perCorrPAC;
                                            par_out(LFPNo).PAC(PACii).which_eventPAC=par_out(LFPrep).PAC(PACii).which_eventPAC;
                                            par_out(LFPNo).PAC(PACii).meanPeakAngle=par_out(LFPrep).PAC(PACii).meanPeakAngle;
                                            par_out(LFPNo).PAC(PACii).meanPeakPower=par_out(LFPrep).PAC(PACii).meanPeakPower;
                                            par_out(LFPNo).PAC(PACii).meanTroughPower=par_out(LFPrep).PAC(PACii).meanTroughPower;
                                            par_out(LFPNo).PAC(PACii).PACtimecourse=par_out(LFPrep).PAC(PACii).PACtimecourse;
                                            
                                            
                                            %Enter the per LFP values for PACwave
                                            par_out(LFPNo).PACwave(PACii).meanPeakPower=par_out(LFPrep).PACwave(PACii).meanPeakPower;
                                            par_out(LFPNo).PACwave(PACii).meanTroughPower=par_out(LFPrep).PACwave(PACii).meanTroughPower;
                                            par_out(LFPNo).PACwave(PACii).meanPower=par_out(LFPrep).PACwave(PACii).meanPower;
                                            par_out(LFPNo).PACwave(PACii).PACtimecourse=par_out(LFPrep).PACwave(PACii).PACtimecourse;
                                            
                                        end
                                    end
                                    fprintf(1, 'For file no %d, LFP no %d the number of trials included in PAC is %d\n',filNum,LFPNo,par_out(LFPNo).PAC(PACii).no_trials);
                                    
                                end
                                
                                
                                fprintf(1, 'WARNING. Number of trials trimmed because For file no %d the number of trials differ between different LFPs\n ',filNum);
                            end
                            
                            t=t_pac;
                            for LFPNo=1:length(handles.drgbchoices.which_electrodes)
                                
                                %                             this_peakLFPNo=handles.drgbchoices.which_electrodes(LFPNo);
                                if LFPNo>1
                                    if last_No_trials~=length(par_out(LFPNo).PACwave(1).meanPeakPower)
                                        fprintf(1, 'WARNING. For file no %d the number of trials differ between different LFPs\n ',filNum);
                                    end
                                end
                                last_No_trials=length(par_out(LFPNo).PACwave(1).meanPeakPower);
                                
                                if LFPNo==1
                                    if no_trialsPACwave==0
                                        all_log_P_timecoursePACwavepeak=zeros(length(handles.drgbchoices.PACburstLowF),length(handles.drgbchoices.which_electrodes),length(par_out(LFPNo).PACwave(1).meanPeakPower),length(t));
                                        all_log_P_timecoursePACwavetrough=zeros(length(handles.drgbchoices.PACburstLowF),length(handles.drgbchoices.which_electrodes),length(par_out(LFPNo).PACwave(1).meanTroughPower),length(t));
                                    else
                                        %Peak
                                        this_all_log_P_timecoursePACwave=[];
                                        this_all_log_P_timecoursePACwave=all_log_P_timecoursePACwavepeak;
                                        szt=size(this_all_log_P_timecoursePACwave,4);
                                        all_log_P_timecoursePACwavepeak=zeros(length(handles.drgbchoices.PACburstLowF),length(handles.drgbchoices.which_electrodes),length(par_out(LFPNo).PACwave(1).meanPeakPower)+no_trialsPACwave,length(t));
                                        all_log_P_timecoursePACwavepeak(:,:,1:no_trialsPACwave,1:szt)=this_all_log_P_timecoursePACwave(:,:,1:no_trialsPACwave,:);
                                        %Trough
                                        this_all_log_P_timecoursePACwave=[];
                                        this_all_log_P_timecoursePACwave=all_log_P_timecoursePACwavetrough;
                                        szt=size(this_all_log_P_timecoursePACwave,4);
                                        all_log_P_timecoursePACwavetrough=zeros(length(handles.drgbchoices.PACburstLowF),length(handles.drgbchoices.which_electrodes),length(par_out(LFPNo).PACwave(1).meanTroughPower)+no_trialsPACwave,length(t));
                                        all_log_P_timecoursePACwavetrough(:,:,1:no_trialsPACwave,1:szt)=this_all_log_P_timecoursePACwave(:,:,1:no_trialsPACwave,:);
                                    end
                                end
                                for PACii=1:length(handles.drgbchoices.PACburstLowF)
                                    for trNo=1:length(par_out(LFPNo).PACwave(1).meanPeakPower)
                                        all_log_P_timecoursePACwavepeak(PACii,LFPNo,no_trialsPACwave+trNo,:)=par_out(LFPNo).PACwave(PACii).PACtimecourse(trNo).peakPower;
                                    end
                                    for trNo=1:length(par_out(LFPNo).PACwave(1).meanTroughPower)
                                        all_log_P_timecoursePACwavetrough(PACii,LFPNo,no_trialsPACwave+trNo,:)=par_out(LFPNo).PACwave(PACii).PACtimecourse(trNo).troughPower;
                                    end
                                end
                                
                                
                            end
                        end
                        
                        %Calculate licks
                        stamped_lick_ii=[];
                        these_stamped_lick_times=[];
                        no_trials_l=[];
                        trials_included_l=[];
                        handles.lastTrialNo=handles.drg.session(1).noTrials;
                        
                        [lick_freq,times_lick_freq,lick_traces,CIlickf,lick_trace_times,stamped_lick_ii...
                            ,these_stamped_lick_times,no_trials_l,trials_included_l, lick_threshold]=drgGetLicks(handles);
                        
                        %Extract the data for powerLFP
                        if (sum(handles.drgbchoices.which_discriminant==1)>0)||(sum(handles.drgbchoices.which_discriminant==2)>0)...
                                ||(sum(handles.drgbchoices.which_discriminant==3)>0)||(sum(handles.drgbchoices.which_discriminant==6)>0)
                            t=par_out(1).t;
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
                        
                        %Extract the data for PAC
                        if (sum(handles.drgbchoices.which_discriminant==4)>0)||(sum(handles.drgbchoices.which_discriminant==5)>0)...
                                ||(sum(handles.drgbchoices.which_discriminant==7)>0)||(sum(handles.drgbchoices.which_discriminant==8)>0)||...
                                (sum(handles.drgbchoices.which_discriminant==9)>0)
                            t=t_power;
                            %Save all_which_events and all_perCorr_pertr
                            if filNum==first_file_for_this_mouse
                                all_which_eventsPAC=zeros(length(handles.drgbchoices.evTypeNos),par_out(1).PAC(1).no_trials);
                                all_which_eventsPAC(:,:)=par_out(1).PAC(1).which_eventPAC;
                                
                                all_perCorr_pertrPAC=zeros(1,par_out(1).PAC(1).no_trials);
                                all_perCorr_pertrPAC(:,:)=par_out(1).PAC(1).perCorrPAC;
                                
                                all_stamped_lick_timesPAC=zeros(par_out(1).PAC(1).no_trials,250);
                                sztslt=size(these_stamped_lick_times);
                                all_stamped_lick_iiPAC=zeros(1,par_out(1).PAC(1).no_trials);
                                for ii=1:par_out(1).PAC(1).no_trials
                                    kk=find(trials_included_l==par_out(1).PAC(1).this_trialNo(ii));
                                    all_stamped_lick_timesPAC(ii,1:sztslt(2))=these_stamped_lick_times(kk,:);
                                    all_stamped_lick_iiPAC(1,ii)=stamped_lick_ii(kk);
                                end
                                
                            else
                                this_all_which_events=[];
                                this_all_which_events=all_which_eventsPAC;
                                all_which_eventsPAC=zeros(length(handles.drgbchoices.evTypeNos),par_out(1).PAC(1).no_trials+no_trialsPAC);
                                all_which_eventsPAC(:,1:no_trialsPAC,:)=this_all_which_events;
                                all_which_eventsPAC(:,no_trialsPAC+1:no_trialsPAC+par_out(1).PAC(1).no_trials,:)=par_out(1).PAC(1).which_eventPAC;
                                
                                this_all_perCorr_pertr=[];
                                this_all_perCorr_pertr=all_perCorr_pertrPAC;
                                all_perCorr_pertrPAC=zeros(1,par_out(1).PAC(1).no_trials+no_trialsPAC);
                                all_perCorr_pertrPAC(:,1:no_trialsPAC)=this_all_perCorr_pertr;
                                all_perCorr_pertrPAC(:,no_trialsPAC+1:no_trialsPAC+par_out(1).PAC(1).no_trials)=par_out(1).PAC(1).perCorrPAC;
                                
                                this_all_stamped_lick_times=[];
                                this_all_stamped_lick_times=all_stamped_lick_timesPAC;
                                all_stamped_lick_timesPAC=zeros(par_out(1).PAC(1).no_trials+no_trialsPAC,250);
                                all_stamped_lick_timesPAC(1:no_trialsPAC,:)=this_all_stamped_lick_times;
                                sztslt=size(these_stamped_lick_times);
                                
                                this_all_stamped_lick_ii=[];
                                this_all_stamped_lick_ii=all_stamped_lick_iiPAC;
                                all_stamped_lick_iiPCA=zeros(1,par_out(1).PAC(1).no_trials+no_trialsPAC);
                                all_stamped_lick_iiPCA(1,1:no_trialsPAC)=this_all_stamped_lick_ii;
                                
                                for ii=1:par_out(1).PAC(1).no_trials
                                    kk=find(trials_included_l==par_out(1).PAC(1).this_trialNo(ii));
                                    all_stamped_lick_times(ii+no_trialsPAC,1:sztslt(2))=these_stamped_lick_times(kk,:);
                                    all_stamped_lick_ii(1,no_trialsPAC+ii)=stamped_lick_ii(1,kk);
                                end
                                
                            end
                            
                            no_trialsPAC=no_trialsPAC+par_out(1).PAC(1).no_trials;
                        end
                        
                        %Extract the data for wavelet power computed at PAC phase
                        if (sum(handles.drgbchoices.which_discriminant==10)>0)||(sum(handles.drgbchoices.which_discriminant==11)>0)||...
                                (sum(handles.drgbchoices.which_discriminant==12)>0)||(sum(handles.drgbchoices.which_discriminant==13)>0)...
                                ||(sum(handles.drgbchoices.which_discriminant==14)>0)||(sum(handles.drgbchoices.which_discriminant==15)>0)...
                                ||(sum(handles.drgbchoices.which_discriminant==16)>0)
                            t=t_pac;
                            
                            if (trials_match==0)||(par_out(1).PAC(1).no_trials~=length(stamped_lick_ii))
                                %Trim the licks
                                
                                these_trials=par_out(1).PAC(1).this_trialNo;
                                old_lick_traces=lick_traces;
                                lick_traces=zeros(length(these_trials),size(lick_traces,2));
                                old_stamped_lick_ii=stamped_lick_ii;
                                stamped_lick_ii=zeros(1,length(these_trials));
                                old_trials_included_l=trials_included_l;
                                trials_included_l=[];
                                
                                for new_trNo=1:par_out(1).PAC(1).no_trials
                                    trNo=find(these_trials(new_trNo)==old_trials_included_l);
                                    lick_traces(new_trNo,:)=old_lick_traces(trNo,:);
                                    stamped_lick_ii(new_trNo)=old_stamped_lick_ii(trNo);
                                    trials_included_l(new_trNo)=old_trials_included_l(trNo);
                                end
                                
                            end
                            
                            %Genetrate the licks
                            these_licks_per_t=zeros(par_out(1).PAC(1).no_trials,length(t));
                            for trNo=1:par_out(1).PAC(1).no_trials
                                for ii_t=1:length(t)
                                    if stamped_lick_ii(trNo)>0
                                        if sum(lick_traces(trNo,(lick_trace_times>(t(ii_t)-((t(2)-t(1))/2)))&(lick_trace_times<=(t(ii_t)+((t(2)-t(1))/2))))>lick_threshold)>0
                                            these_licks_per_t(trNo,ii_t)=1;
                                        end
                                    end
                                end
                            end
                            
                            %Save all_which_events and all_perCorr_pertr
                            if filNum==first_file_for_this_mouse
                                all_which_eventsPACwave=zeros(length(handles.drgbchoices.evTypeNos),par_out(1).PAC(1).no_trials);
                                all_which_eventsPACwave(:,:)=par_out(1).PAC(1).which_eventPAC;
                                
                                all_perCorr_pertrPACwave=zeros(1,par_out(1).PAC(1).no_trials);
                                all_perCorr_pertrPACwave(:,:)=par_out(1).PAC(1).perCorrPAC;
                                
                                all_groups_pertrPACwave=groupNo*ones(1,par_out(1).PAC(1).no_trials);
                                
                                %                                 all_stamped_lick_timesPACwave=zeros(par_out(1).PAC(1).no_trials,250);
                                %                                 sztslt=size(these_stamped_lick_times);
                                all_stamped_lick_iiPACwave=zeros(1,par_out(1).PAC(1).no_trials);
                                all_stamped_lick_iiPACwave=stamped_lick_ii;
                                
                                %                                 for ii=1:par_out(1).PAC(1).no_trials
                                %                                     kk=find(trials_included_l==par_out(1).PAC(1).this_trialNo(ii));
                                %                                     all_stamped_lick_timesPACwave(ii,1:sztslt(2))=these_stamped_lick_times(kk,:);
                                %                                     all_stamped_lick_iiPACwave(1,ii)=stamped_lick_ii(kk);
                                %                                 end
                                
                                all_licks_per_tPACwave=zeros(par_out(1).PAC(1).no_trials,length(t));
                                all_licks_per_tPACwave(:,:)=these_licks_per_t;
                                
                            else
                                this_all_which_eventswave=[];
                                this_all_which_eventswave=all_which_eventsPACwave;
                                all_which_eventsPACwave=zeros(length(handles.drgbchoices.evTypeNos),par_out(1).PAC(1).no_trials+no_trialsPACwave);
                                all_which_eventsPACwave(:,1:no_trialsPACwave,:)=this_all_which_eventswave;
                                all_which_eventsPACwave(:,no_trialsPACwave+1:no_trialsPACwave+par_out(1).PAC(1).no_trials,:)=par_out(1).PAC(1).which_eventPAC;
                                
                                this_all_perCorr_pertrwave=[];
                                this_all_perCorr_pertrwave=all_perCorr_pertrPACwave;
                                all_perCorr_pertrPACwave=zeros(1,par_out(1).PAC(1).no_trials+no_trialsPACwave);
                                all_perCorr_pertrPACwave(:,1:no_trialsPACwave)=this_all_perCorr_pertrwave;
                                all_perCorr_pertrPACwave(:,no_trialsPACwave+1:no_trialsPACwave+par_out(1).PAC(1).no_trials)=par_out(1).PAC(1).perCorrPAC;
                                
                                this_all_groups_pertrwave=[];
                                this_all_groups_pertrwave=all_groups_pertrPACwave;
                                all_groups_pertrPACwave=zeros(1,par_out(1).PAC(1).no_trials+no_trialsPACwave);
                                all_groups_pertrPACwave(:,1:no_trialsPACwave)=this_all_groups_pertrwave;
                                all_groups_pertrPACwave(:,no_trialsPACwave+1:no_trialsPACwave+par_out(1).PAC(1).no_trials)=groupNo*ones(1,par_out(1).PAC(1).no_trials);
                                
                                this_all_stamped_lick_iiPACwave=[];
                                this_all_stamped_lick_iiPACwave=all_stamped_lick_iiPACwave;
                                all_stamped_lick_iiPACwave=zeros(1,par_out(1).PAC(1).no_trials+no_trialsPACwave);
                                all_stamped_lick_iiPACwave(1,1:no_trialsPACwave)=this_all_stamped_lick_iiPACwave;
                                all_stamped_lick_iiPACwave(1,no_trialsPACwave+1:no_trialsPACwave+par_out(1).PAC(1).no_trials)=stamped_lick_ii;
                                
                                this_all_licks_per_t=[];
                                this_all_licks_per_t=all_licks_per_tPACwave;
                                all_licks_per_tPACwave=zeros(par_out(1).PAC(1).no_trials+no_trialsPACwave,length(t));
                                all_licks_per_tPACwave(1:no_trialsPACwave,:)=this_all_licks_per_t;
                                all_licks_per_tPACwave(no_trialsPACwave+1:no_trialsPACwave+par_out(1).PAC(1).no_trials,:)=these_licks_per_t;
                                
                            end
                            
                            no_trialsPACwave=no_trialsPACwave+par_out(1).PAC(1).no_trials;
                        end
                        
                        
                        %Delete temp files
                        
                        
                        %Delete temp files
                        delete([jtPathName jtFileName]);
                        delete([jtPathName jtFileName(10:end-4) '_drg.mat'])
                        
                        if exist([jtPathNameSource this_rhd])==0
                            this_dg=[jtFileName(10:end-4) '.dg'];
                            delete([jtPathName this_dg]);
                            %                                 handles.drg.drta_p.fullName=[jtPathName jtFileName(10:end-4) '.dg'];
                        else
                            this_rhd=[jtFileName(10:end-4) '.rhd'];
                            delete([jtPathName this_rhd]);
                            %                                 handles.drg.drta_p.fullName=[jtPathName jtFileName(10:end-4) '.rhd'];
                        end
                    end
                end
                
                
                %Do discriminant analysis for the percent correct windows
                if mouse_has_data==1
                    
                    %Calculate discriminant analysis and PCA for LFP power
                    if (sum(handles.drgbchoices.which_discriminant==1)>0)||(sum(handles.drgbchoices.which_discriminant==2)>0)...
                            ||(sum(handles.drgbchoices.which_discriminant==3)>0)||(sum(handles.drgbchoices.which_discriminant==6)>0)
                        t=par_out(1).t;
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
                                        %Linear discriminant analysis for
                                        %LFP power
                                        
                                        %Number of trials
                                        %                                         N=sum(these_per_corr);
                                        
                                        these_all_log_P_timecourse=zeros(length(handles.drgbchoices.which_electrodes),sum(these_per_corr),length(t));
                                        these_all_log_P_timecourse(:,:,:)=all_log_P_timecourse(bwii,:,these_per_corr,:);
                                        
                                        these_all_which_events=zeros(length(handles.drgbchoices.events_to_discriminate),sum(these_per_corr));
                                        for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                            kk=find(handles.drgbchoices.evTypeNos==handles.drgbchoices.events_to_discriminate(ii));
                                            these_all_which_events(ii,:)= all_which_events(kk,these_per_corr);
                                        end
                                        
                                        
                                        
                                        fprintf(1, ['LDA processed for LFP power for mouse No %d ' handles.drgbchoices.group_no_names{groupNo} ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' with %d trials\n'],mouseNo,N);
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
                                        dimensionality=zeros(1,length(t));
                                        auROC=zeros(1,length(t));
                                        per_targets=zeros(length(handles.drgbchoices.events_to_discriminate),N);
                                        
                                        parfor time_point=1:length(t)
                                            
                                            %LFP power per trial per electrode
                                            measurements=zeros(N,length(handles.drgbchoices.which_electrodes));
                                            measurements(:,:)=these_all_log_P_timecourse(:,:,time_point)';
                                            
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
                                    
                                    
                                    if sum(handles.drgbchoices.which_discriminant==6)>0
                                        %Linear discriminant analysis for
                                        %subsets of electrodes for LFP
                                        %power
                                        t_from=2.15;
                                        t_to=2.5;
                                        subt=t((t>=t_from)&(t<=t_to));
                                        
                                        these_all_log_P_timecourse=zeros(length(handles.drgbchoices.which_electrodes),sum(these_per_corr),length(subt));
                                        these_all_log_P_timecourse(:,:,:)=all_log_P_timecourse(bwii,:,these_per_corr,(t>=t_from)&(t<=t_to));
                                        
                                        these_all_which_events=zeros(length(handles.drgbchoices.events_to_discriminate),sum(these_per_corr));
                                        for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                            kk=find(handles.drgbchoices.evTypeNos==handles.drgbchoices.events_to_discriminate(ii));
                                            these_all_which_events(ii,:)= all_which_events(kk,these_per_corr);
                                        end
                                        
                                        fprintf(1, ['LDA processed for LFP power for mouse No %d ' handles.drgbchoices.group_no_names{groupNo} ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' with %d trials\n'],mouseNo,N);
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
                                                    measurements(:,:)=these_all_log_P_timecourse(elect_combs(noelc,:),:,time_point)';
                                                    
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
                                                    fprintf(1, 'LDA percent correct classification %d (for timepoint %d and number of electrodes %d)\n',100*sum(per_targets(1,:)==test_out(1,:))/N,time_point,no_elect);
                                                    
                                                end
                                            end
                                        end
                                        
                                        if bwii==4
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
                                            
                                            suptitle(['LFP power LDA ' handles.drgbchoices.bwlabels{bwii} ' mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                                        end
                                        
                                        handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).ecomb_discriminant_calculated=1;
                                        handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).ecomb_par_out=par_out;
                                    end
                                    
                                    if (sum(handles.drgbchoices.which_discriminant==1)>0)||(sum(handles.drgbchoices.which_discriminant==2)>0)
                                        
                                        
                                        if bwii==4
                                            figNo=figNo+1;
                                            try
                                                close(figNo)
                                            catch
                                            end
                                            
                                            figure(figNo)
                                            
                                            subplot(1,3,1)
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
                                            
                                            subplot(1,3,2)
                                            hold on
                                            
                                            plot(t,auROC)
                                            
                                            %Odor on markers
                                            plot([0 0],[-0.3 0.5],'-k')
                                            odorhl=plot([0 2.5],[0 0],'-k','LineWidth',5);
                                            plot([2.5 2.5],[-0.3 0.5],'-k')
                                            
                                            %title(['auROC for LDA for ' handles.drgbchoices.bwlabels{bwii} ' mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                                            xlabel('Time (sec)')
                                            ylabel('auROC')
                                            ylim([-0.3 0.6])
                                            
                                            %Dimensionality
                                            subplot(1,3,3)
                                            hold on
                                            
                                            plot(t,dimensionality)
                                            
                                            %Odor on markers
                                            maxdim=max(dimensionality);
                                            mindim=min(dimensionality);
                                            plot([0 0],[mindim-0.1*(maxdim-mindim) maxdim+0.1*(maxdim-mindim)],'-k')
                                            odorhl=plot([0 2.5],[mindim mindim],'-k','LineWidth',5);
                                            plot([2.5 2.5],[mindim-0.1*(maxdim-mindim) maxdim+0.1*(maxdim-mindim)],'-k')
                                            
                                            %title(['auROC for LDA for ' handles.drgbchoices.bwlabels{bwii} ' mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                                            xlabel('Time (sec)')
                                            ylabel('Dimensionality')
                                            
                                            suptitle(['LFP power LDA analysis for ' handles.drgbchoices.bwlabels{bwii} ' mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                                        end
                                        
                                        handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated=1;
                                        handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).discriminant_correct=zeros(1,length(t));
                                        handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).discriminant_correct(1,:)=discriminant_correct(1,:);
                                        
                                        handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).dimensionality=zeros(1,length(t));
                                        handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).dimensionality(1,:)=dimensionality(1,:);
                                        
                                        
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
                                        
                                        
                                        
                                        fprintf(1, ['PCA processed for LFP power for mouse No %d ' handles.drgbchoices.group_no_names{groupNo} ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' with %d trials\n'],mouseNo,N);
                                        fprintf(1,'For these events: ')
                                        for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                            fprintf(1,[handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(ii)} ' '])
                                        end
                                        fprintf(1,'\n')
                                        
                                        for time_point=1:length(t)
                                            par_t_out(time_point).principal_components=zeros(N,length(handles.drgbchoices.which_electrodes));
                                        end
                                        
                                        
                                        
                                        %                                         test_out_per_timepoint=zeros(length(handles.drgbchoices.events_to_discriminate),N,length(t));
                                        %                                         shuffled_out_per_timepoint=zeros(length(handles.drgbchoices.events_to_discriminate),N,length(t));
                                        %
                                        for time_point=1:length(t)
                                            par_t_out(time_point).principal_components=[];
                                        end
                                        
                                        for time_point=1:length(t)
                                            
                                            %LFP power per trial per electrode
                                            measurements=zeros(N,length(handles.drgbchoices.which_electrodes));
                                            measurements(:,:)=these_all_log_P_timecourse(:,:,time_point)';
                                            
                                            %Do the PCA
                                            [coeff,par_t_out(time_point).principal_components,par_t_out(time_point).PC_variance]=pca(measurements);
                                            
                                        end
                                        
                                        
                                        %Show a figure of the PCA and record the output
                                        
                                        principal_components=zeros(length(t),N,length(handles.drgbchoices.which_electrodes));
                                        PC_variance=zeros(length(t),length(handles.drgbchoices.which_electrodes));
                                        for time_point=1:length(t)
                                            principal_components(time_point,:,:)=par_t_out(time_point).principal_components;
                                            PC_variance(time_point,:)=par_t_out(time_point).PC_variance;
                                        end
                                        
                                        if bwii==4
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
                                            title(['LFP power PC1 for ' handles.drgbchoices.bwlabels{bwii} ' mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                                            xlabel('Time (sec)')
                                            ylabel('PC1')
                                        end
                                        
                                        handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated=1;
                                        handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).principal_components=zeros(length(t),N,length(handles.drgbchoices.which_electrodes));
                                        handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).principal_components=principal_components;
                                        handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).PC_variance=zeros(length(t),length(handles.drgbchoices.which_electrodes));
                                        handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).PC_variance=PC_variance;
                                        
                                        
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
                    end
                    
                    
                    %Calculate discriminant analysis for hogh frequency envelope power calculated at different phases
                    if (sum(handles.drgbchoices.which_discriminant==7)>0)||(sum(handles.drgbchoices.which_discriminant==8)>0)...
                            ||(sum(handles.drgbchoices.which_discriminant==9)>0)
                        t=t_power;
                        for PACii=1:length(handles.drgbchoices.PACburstLowF)
                            
                            for percent_correct_ii=1:length(handles.drgbchoices.per_lab)
                                
                                discriminant_correct=zeros(1,length(t));
                                discriminant_correct_shuffled=zeros(1,length(t));
                                auROC=zeros(1,length(t));
                                
                                these_per_corr=(all_perCorr_pertrPAC>=handles.drgbchoices.percent_windows(percent_correct_ii,1))...
                                    &(all_perCorr_pertrPAC<=handles.drgbchoices.percent_windows(percent_correct_ii,2));
                                N=sum(these_per_corr);
                                
                                %Do the analysis only if there are more than 20 trials
                                if N>=20
                                    
                                    if (sum(handles.drgbchoices.which_discriminant==7)>0)
                                        %Linear discriminant analysis for
                                        %peak PAC power
                                        
                                        these_all_log_P_timecoursePAC=zeros(length(handles.drgbchoices.which_electrodes),sum(these_per_corr),length(t));
                                        these_all_log_P_timecoursePAC(:,:,:)=all_log_P_timecoursePACpeak(PACii,:,these_per_corr,:);
                                        
                                        these_all_which_events=zeros(length(handles.drgbchoices.events_to_discriminate),sum(these_per_corr));
                                        for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                            kk=find(handles.drgbchoices.evTypeNos==handles.drgbchoices.events_to_discriminate(ii));
                                            these_all_which_events(ii,:)= all_which_eventsPAC(kk,these_per_corr);
                                        end
                                        
                                        
                                        
                                        fprintf(1, ['LDA processed for peak PAC power for mouse No %d ' handles.drgbchoices.group_no_names{groupNo} ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' with %d trials\n'],mouseNo,N);
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
                                            
                                            %LFP power per trial per electrode
                                            measurements=zeros(N,length(handles.drgbchoices.which_electrodes));
                                            measurements(:,:)=these_all_log_P_timecoursePAC(:,:,time_point)';
                                            
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
                                            fprintf(1, 'LDA for peak PAC power percent correct classification %d (for timepoint %d out of %d)\n',100*sum(per_targets(1,:)==test_out(1,:))/N,time_point,length(t));
                                            
                                        end
                                        
                                        
                                        if PACii==3
                                            figNo=figNo+1;
                                            try
                                                close(figNo)
                                            catch
                                            end
                                            
                                            figure(figNo)
                                            
                                            subplot(2,3,1)
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
                                            
                                            subplot(2,3,2)
                                            hold on
                                            
                                            plot(t,auROC)
                                            
                                            %Odor on markers
                                            plot([0 0],[-0.3 0.5],'-k')
                                            odorhl=plot([0 2.5],[0 0],'-k','LineWidth',5);
                                            plot([2.5 2.5],[0 0.5],'-k')
                                            
                                            %title(['auROC for LDA for ' handles.drgbchoices.bwlabels{PACii} ' mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                                            xlabel('Time (sec)')
                                            ylabel('auROC peak')
                                            ylim([-0.3 0.6])
                                            
                                            
                                            subplot(2,3,3)
                                            hold on
                                            
                                            plot(t,dimensionality)
                                            
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
                                        end
                                        %suptitle(['PAC power LDA analysis for Theta/' handles.drgbchoices.PACnames{PACii} ' PAC, mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                                        
                                        handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated=1;
                                        handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_peak=zeros(1,length(t));
                                        handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_peak(1,:)=discriminant_correct(1,:);
                                        
                                        handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).dimensionality_peak=zeros(1,length(t));
                                        handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).dimensionality_peak(1,:)=dimensionality(1,:);
                                        
                                        
                                        handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).auROC_peak=zeros(1,length(t));
                                        handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).auROC_peak=auROC;
                                        
                                        handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_shuffled_peak=zeros(1,length(t));
                                        handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_shuffled_peak(1,:)=discriminant_correct_shuffled(1,:);
                                        
                                        handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).test_out_per_timepoint_peak=test_out_per_timepoint;
                                        handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).shuffled_out_per_timepoint_peak=shuffled_out_per_timepoint;
                                        
                                        handles_out.t_power=t';
                                        handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).no_trials=N;
                                        handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).per_targets=per_targets;
                                        
                                        these_all_which_events=[];
                                        these_all_which_events=all_which_eventsPAC(:,these_per_corr);
                                        handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).which_events=these_all_which_events;
                                        
                                        %                                         these_all_stamped_lick_times=[];
                                        %                                         these_all_stamped_lick_times=all_stamped_lick_times(these_per_corr,:);
                                        %                                         handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).stamped_lick_times=these_all_stamped_lick_times;
                                        %
                                        %                                         these_all_stamped_lick_ii=[];
                                        %                                         these_all_stamped_lick_ii=all_stamped_lick_ii(1,these_per_corr);
                                        %                                         handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).stamped_lick_ii=these_all_stamped_lick_ii;
                                        
                                        %Linear discriminant analysis for
                                        %trough PAC power
                                        
                                        these_all_log_P_timecoursePAC=zeros(length(handles.drgbchoices.which_electrodes),sum(these_per_corr),length(t));
                                        these_all_log_P_timecoursePAC(:,:,:)=all_log_P_timecoursePACtrough(PACii,:,these_per_corr,:);
                                        
                                        these_all_which_events=zeros(length(handles.drgbchoices.events_to_discriminate),sum(these_per_corr));
                                        for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                            kk=find(handles.drgbchoices.evTypeNos==handles.drgbchoices.events_to_discriminate(ii));
                                            these_all_which_events(ii,:)= all_which_eventsPAC(kk,these_per_corr);
                                        end
                                        
                                        
                                        
                                        fprintf(1, ['LDA processed for trough PAC power for mouse No %d ' handles.drgbchoices.group_no_names{groupNo} ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' with %d trials\n'],mouseNo,N);
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
                                            
                                            %LFP power per trial per electrode
                                            measurements=zeros(N,length(handles.drgbchoices.which_electrodes));
                                            measurements(:,:)=these_all_log_P_timecoursePAC(:,:,time_point)';
                                            
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
                                            fprintf(1, 'LDA for trough PAC power percent correct classification %d (for timepoint %d out of %d)\n',100*sum(per_targets(1,:)==test_out(1,:))/N,time_point,length(t));
                                            
                                        end
                                        
                                        if PACii==3
                                            subplot(2,3,4)
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
                                            ylabel('% correct trough')
                                            
                                            subplot(2,3,5)
                                            hold on
                                            
                                            plot(t,auROC)
                                            
                                            %Odor on markers
                                            plot([0 0],[0 0.5],'-k')
                                            odorhl=plot([0 2.5],[0 0],'-k','LineWidth',5);
                                            plot([2.5 2.5],[0 0.5],'-k')
                                            
                                            %title(['auROC for LDA for ' handles.drgbchoices.bwlabels{PACii} ' mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                                            xlabel('Time (sec)')
                                            ylabel('auROC trough')
                                            ylim([-0.3 0.6])
                                            
                                            subplot(2,3,6)
                                            hold on
                                            
                                            plot(t,dimensionality)
                                            
                                            mindim=min(dimensionality);
                                            maxdim=max(dimensionality);
                                            
                                            %Odor on markers
                                            plot([0 0],[mindim-0.1*(maxdim-mindim) maxdim+0.1*(maxdim-mindim)],'-k')
                                            odorhl=plot([0 2.5],[0 0],'-k','LineWidth',5);
                                            plot([2.5 2.5],[mindim-0.1*(maxdim-mindim) maxdim+0.1*(maxdim-mindim)],'-k')
                                            
                                            %title(['auROC for LDA for ' handles.drgbchoices.bwlabels{PACii} ' mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                                            xlabel('Time (sec)')
                                            ylabel('Dimensionality trough')
                                            
                                            suptitle(['PAC power LDA analysis for Theta/' handles.drgbchoices.PACnames{PACii} ' PAC, mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                                        end
                                        
                                        handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated=1;
                                        handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_trough=zeros(1,length(t));
                                        handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_trough(1,:)=discriminant_correct(1,:);
                                        
                                        handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).dimensionality_trough=zeros(1,length(t));
                                        handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).dimensionality_trough(1,:)=dimensionality(1,:);
                                        
                                        handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).auROC_trough=zeros(1,length(t));
                                        handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).auROC_trough=auROC;
                                        
                                        handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_shuffled_trough=zeros(1,length(t));
                                        handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_shuffled_trough(1,:)=discriminant_correct_shuffled(1,:);
                                        
                                        handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).test_out_per_timepoint_trough=test_out_per_timepoint;
                                        handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).shuffled_out_per_timepoint_trough=shuffled_out_per_timepoint;
                                        
                                        %                                         handles_out.t_power=t';
                                        %                                         handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).no_trials=N;
                                        %                                         handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).per_targets=per_targets;
                                        %
                                        %                                         these_all_which_events=[];
                                        %                                         these_all_which_events=all_which_eventsPAC(:,these_per_corr);
                                        %                                         handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).which_events=these_all_which_events;
                                        %
                                        %                                         these_all_stamped_lick_times=[];
                                        %                                         these_all_stamped_lick_times=all_stamped_lick_times(these_per_corr,:);
                                        %                                         handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).stamped_lick_times=these_all_stamped_lick_times;
                                        %
                                        %                                         these_all_stamped_lick_ii=[];
                                        %                                         these_all_stamped_lick_ii=all_stamped_lick_ii(1,these_per_corr);
                                        %                                         handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).stamped_lick_ii=these_all_stamped_lick_ii;
                                        
                                        
                                    end
                                    
                                    if sum(handles.drgbchoices.which_discriminant==8)>0
                                        
                                        %PCA for peak PAC power
                                        
                                        these_all_log_P_timecoursePAC=zeros(length(handles.drgbchoices.which_electrodes),sum(these_per_corr),length(t));
                                        these_all_log_P_timecoursePAC(:,:,:)=all_log_P_timecoursePACpeak(PACii,:,these_per_corr,:);
                                        
                                        these_all_which_events=zeros(length(handles.drgbchoices.events_to_discriminate),sum(these_per_corr));
                                        for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                            kk=find(handles.drgbchoices.evTypeNos==handles.drgbchoices.events_to_discriminate(ii));
                                            these_all_which_events(ii,:)= all_which_eventsPAC(kk,these_per_corr);
                                        end
                                        
                                        fprintf(1, ['PCA processed for peak PAC power for mouse No %d ' handles.drgbchoices.group_no_names{groupNo} ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' with %d trials\n'],mouseNo,N);
                                        fprintf(1,'For these events: ')
                                        for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                            fprintf(1,[handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(ii)} ' '])
                                        end
                                        fprintf(1,'\n')
                                        
                                        for time_point=1:length(t)
                                            par_t_out(time_point).principal_components=zeros(N,length(handles.drgbchoices.which_electrodes));
                                        end
                                        
                                        
                                        for time_point=1:length(t)
                                            
                                            %LFP power per trial per electrode
                                            measurements=zeros(N,length(handles.drgbchoices.which_electrodes));
                                            measurements(:,:)=these_all_log_P_timecoursePAC(:,:,time_point)';
                                            
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
                                        
                                        principal_components=zeros(length(t),N,length(handles.drgbchoices.which_electrodes));
                                        PC_variance=zeros(length(t),length(handles.drgbchoices.which_electrodes));
                                        for time_point=1:length(t)
                                            principal_components(time_point,:,:)=par_t_out(time_point).principal_components;
                                            PC_variance(time_point,:)=par_t_out(time_point).PC_variance;
                                        end
                                        
                                        %Show the result of the PCA
                                        if PACii==3
                                            figNo=figNo+1;
                                            try
                                                close(figNo)
                                            catch
                                            end
                                            
                                            figure(figNo)
                                            
                                            %Show PCA before odor on
                                            subplot(2,4,5)
                                            hold on
                                            
                                            these_pcs=zeros(N,length(handles.drgbchoices.which_electrodes));
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
                                            
                                            these_pcs=zeros(N,length(handles.drgbchoices.which_electrodes));
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
                                        end
                                        
                                        handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PCA_calculated=1;
                                        handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).principal_components_peak=zeros(length(t),N,length(handles.drgbchoices.which_electrodes));
                                        handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).principal_components_peak=principal_components;
                                        handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).PC_variance_peak=zeros(length(t),length(handles.drgbchoices.which_electrodes));
                                        handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).PC_variance_peak=PC_variance;
                                        
                                        
                                        handles_out.t_power=t_power';
                                        handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).no_trials=N;
                                        
                                        these_all_which_events=[];
                                        these_all_which_events=all_which_eventsPAC(:,these_per_corr);
                                        handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).which_events=these_all_which_events;
                                        
                                        %                                         these_all_stamped_lick_times=[];
                                        %                                         these_all_stamped_lick_times=all_stamped_lick_times(these_per_corr,:);
                                        %                                         handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).stamped_lick_times=these_all_stamped_lick_times;
                                        %
                                        %                                         these_all_stamped_lick_ii=[];
                                        %                                         these_all_stamped_lick_ii=all_stamped_lick_ii(1,these_per_corr);
                                        %                                         handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).stamped_lick_ii=these_all_stamped_lick_ii;
                                        
                                        
                                        %PCA for trough PAC power
                                        
                                        these_all_log_P_timecoursePAC=zeros(length(handles.drgbchoices.which_electrodes),sum(these_per_corr),length(t));
                                        these_all_log_P_timecoursePAC(:,:,:)=all_log_P_timecoursePACtrough(PACii,:,these_per_corr,:);
                                        
                                        these_all_which_events=zeros(length(handles.drgbchoices.events_to_discriminate),sum(these_per_corr));
                                        for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                            kk=find(handles.drgbchoices.evTypeNos==handles.drgbchoices.events_to_discriminate(ii));
                                            these_all_which_events(ii,:)= all_which_eventsPAC(kk,these_per_corr);
                                        end
                                        
                                        fprintf(1, ['PCA processed for trough PAC power for mouse No %d ' handles.drgbchoices.group_no_names{groupNo} ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' with %d trials\n'],mouseNo,N);
                                        fprintf(1,'For these events: ')
                                        for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                            fprintf(1,[handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(ii)} ' '])
                                        end
                                        fprintf(1,'\n')
                                        
                                        for time_point=1:length(t)
                                            par_t_out(time_point).principal_components=zeros(N,length(handles.drgbchoices.which_electrodes));
                                        end
                                        
                                        for time_point=1:length(t)
                                            par_t_out(time_point).principal_components=[];
                                        end
                                        
                                        for time_point=1:length(t)
                                            
                                            %LFP power per trial per electrode
                                            measurements=zeros(N,length(handles.drgbchoices.which_electrodes));
                                            measurements(:,:)=these_all_log_P_timecoursePAC(:,:,time_point)';
                                            
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
                                        
                                        principal_components=zeros(length(t),N,length(handles.drgbchoices.which_electrodes));
                                        PC_variance=zeros(length(t),length(handles.drgbchoices.which_electrodes));
                                        for time_point=1:length(t)
                                            principal_components(time_point,:,:)=par_t_out(time_point).principal_components;
                                            PC_variance(time_point,:)=par_t_out(time_point).PC_variance;
                                        end
                                        
                                        
                                        if PACii==3
                                            %Show PCA before odor on
                                            subplot(2,4,7)
                                            hold on
                                            
                                            these_pcs=zeros(N,length(handles.drgbchoices.which_electrodes));
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
                                            
                                            these_pcs=zeros(N,length(handles.drgbchoices.which_electrodes));
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
                                            suptitle(['PAC power PC1 for Theta/' handles.drgbchoices.PACnames{PACii} ' PAC, for mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                                            title('PC1 for trough')
                                            xlabel('Time (sec)')
                                            ylabel('PC1')
                                        end
                                        
                                        %                                         handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated=1;
                                        handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).principal_components_trough=zeros(length(t),N,length(handles.drgbchoices.which_electrodes));
                                        handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).principal_components_trough=principal_components;
                                        handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).PC_variance_trough=zeros(length(t),length(handles.drgbchoices.which_electrodes));
                                        handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).PC_variance_trough=PC_variance;
                                        
                                        %                                         handles_out.t_power=t';
                                        %                                         handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).no_trials=N;
                                        %
                                        %                                         these_all_which_events=[];
                                        %                                         these_all_which_events=all_which_eventsPAC(:,these_per_corr);
                                        %                                         handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).which_events=these_all_which_events;
                                        %
                                        %                                         these_all_stamped_lick_times=[];
                                        %                                         these_all_stamped_lick_times=all_stamped_lick_times(these_per_corr,:);
                                        %                                         handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).stamped_lick_times=these_all_stamped_lick_times;
                                        %
                                        %                                         these_all_stamped_lick_ii=[];
                                        %                                         these_all_stamped_lick_ii=all_stamped_lick_ii(1,these_per_corr);
                                        %                                         handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).stamped_lick_ii=these_all_stamped_lick_ii;
                                        
                                    end
                                    
                                    if sum(handles.drgbchoices.which_discriminant==9)>0
                                        %Linear discriminant analysis for
                                        %subsets of electrodes for LFP
                                        %power
                                        t_from=2.15;
                                        t_to=2.5;
                                        subt=t((t>=t_from)&(t<=t_to));
                                        
                                        these_all_log_P_timecoursePAC=zeros(length(handles.drgbchoices.which_electrodes),sum(these_per_corr),length(subt));
                                        these_all_log_P_timecoursePAC(:,:,:)=all_log_P_timecoursePAC(PACii,:,these_per_corr,(t>=t_from)&(t<=t_to));
                                        
                                        these_all_which_events=zeros(length(handles.drgbchoices.events_to_discriminate),sum(these_per_corr));
                                        for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                            kk=find(handles.drgbchoices.evTypeNos==handles.drgbchoices.events_to_discriminate(ii));
                                            these_all_which_events(ii,:)= all_which_eventsPAC(kk,these_per_corr);
                                        end
                                        
                                        fprintf(1, ['LDA processed for PAC power for mouse No %d ' handles.drgbchoices.group_no_names{groupNo} ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' with %d trials\n'],mouseNo,N);
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
                                                    measurements(:,:)=these_all_log_P_timecoursePAC(elect_combs(noelc,:),:,time_point)';
                                                    
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
                                                    fprintf(1, 'LDA percent correct classification %d (for timepoint %d and number of electrodes %d)\n',100*sum(per_targets(1,:)==test_out(1,:))/N,time_point,no_elect);
                                                    
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
                                            
                                            suptitle(['PAC power LDA ' handles.drgbchoices.bwlabels{PACii} ' mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                                        end
                                        
                                        handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).ecomb_discriminant_calculated=1;
                                        handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).ecomb_par_out=par_out;
                                    end
                                    
                                else
                                    handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated=0;
                                    fprintf(1, ['LDA/PCA not processed for PAC power for mouse No %d ' handles.drgbchoices.group_no_names{groupNo} ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' because there were only %d trials (fewer than 20 trials)\n'],mouseNo,N);
                                end
                            end
                        end
                    end
                    
                    %Calculate discriminant analysis for wavelet power calculated at different PAC phases
                    if (sum(handles.drgbchoices.which_discriminant==10)>0)||(sum(handles.drgbchoices.which_discriminant==11)>0)...
                            ||(sum(handles.drgbchoices.which_discriminant==12)>0)
                        t=t_pac;
                        for PACii=1:length(handles.drgbchoices.PACburstLowF)
                            
                            for percent_correct_ii=1:length(handles.drgbchoices.per_lab)
                                
                                discriminant_correct=zeros(1,length(t));
                                discriminant_correct_shuffled=zeros(1,length(t));
                                auROC=zeros(1,length(t));
                                
                                these_per_corr=(all_perCorr_pertrPACwave>=handles.drgbchoices.percent_windows(percent_correct_ii,1))...
                                    &(all_perCorr_pertrPACwave<=handles.drgbchoices.percent_windows(percent_correct_ii,2));
                                N=sum(these_per_corr);
                                
                                
                                %Stamp the licks
                                these_all_licks_per_tPACwave=[];
                                these_all_licks_per_tPACwave=all_licks_per_tPACwave(these_per_corr,:);
                                
                                these_all_stamped_lick_ii=[];
                                these_all_stamped_lick_ii=all_stamped_lick_iiPACwave(1,these_per_corr);
                                
                                
                                %Do the analysis only if there are more than 20 trials
                                if N>=20
                                    
                                    if (sum(handles.drgbchoices.which_discriminant==10)>0)
                                        %Linear discriminant analysis for
                                        %peak PAC power
                                        
                                        these_all_log_P_timecoursePACwave=zeros(length(handles.drgbchoices.which_electrodes),sum(these_per_corr),length(t));
                                        these_all_log_P_timecoursePACwave(:,:,:)=all_log_P_timecoursePACwavepeak(PACii,:,these_per_corr,:);
                                        
                                        these_all_which_events=zeros(length(handles.drgbchoices.events_to_discriminate),sum(these_per_corr));
                                        for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                            kk=find(handles.drgbchoices.evTypeNos==handles.drgbchoices.events_to_discriminate(ii));
                                            these_all_which_events(ii,:)= all_which_eventsPACwave(kk,these_per_corr);
                                        end
                                        
                                        
                                        
                                        fprintf(1, ['LDA processed for peak PAC wavelet power for mouse No %d ' handles.drgbchoices.group_no_names{groupNo} ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' with %d trials\n'],mouseNo,N);
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
                                            
                                            %LFP power per trial per electrode
                                            measurements=zeros(N,length(handles.drgbchoices.which_electrodes));
                                            measurements(:,:)=these_all_log_P_timecoursePACwave(:,:,time_point)';
                                            
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
                                            fprintf(1, 'LDA for peak PAC wavelet power percent correct classification %d (for timepoint %d out of %d)\n',100*sum(per_targets(1,:)==test_out(1,:))/N,time_point,length(t));
                                            
                                        end
                                        
                                        
                                        if PACii==3
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
                                        end
                                        
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
                                            splus_licks=zeros(1,sum((these_all_which_events(1,:)==1)&(these_all_stamped_lick_ii>0)));
                                            splus_licks(1,:)=these_all_licks_per_tPACwave((these_all_which_events(1,:)==1)&(these_all_stamped_lick_ii>0),ii_t);
                                            sminus_licks=zeros(1,sum((these_all_which_events(1,:)==0)&(these_all_stamped_lick_ii>0)));
                                            sminus_licks(1,:)=these_all_licks_per_tPACwave((these_all_which_events(1,:)==0)&(these_all_stamped_lick_ii>0),ii_t);
                                            p_val_lick(ii_t)=ranksum(splus_licks,sminus_licks);
                                        end
                                        
                                        
                                        %suptitle(['PAC power LDA analysis for Theta/' handles.drgbchoices.PACnames{PACii} ' PAC, mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).all_stamped_lick_ii=these_all_stamped_lick_ii;
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).all_licks_per_tPACwave=these_all_licks_per_tPACwave;
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).p_val_peak=p_val_peak;
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).p_val_lick=p_val_lick;
                                        
                                        handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated=1;
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_peak=zeros(1,length(t));
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_peak(1,:)=discriminant_correct(1,:);
                                        
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).dimensionality_peak=zeros(1,length(t));
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).dimensionality_peak(1,:)=dimensionality(1,:);
                                        
                                        
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).auROC_peak=zeros(1,length(t));
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).auROC_peak=auROC;
                                        
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_shuffled_peak=zeros(1,length(t));
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_shuffled_peak(1,:)=discriminant_correct_shuffled(1,:);
                                        
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).test_out_per_timepoint_peak=test_out_per_timepoint;
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).shuffled_out_per_timepoint_peak=shuffled_out_per_timepoint;
                                        
                                        handles_out.t_power=t';
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).no_trials=N;
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).per_targets=per_targets;
                                        
                                        these_all_which_events=[];
                                        these_all_which_events=all_which_eventsPACwave(:,these_per_corr);
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).which_events=these_all_which_events;
                                        
                                        %                                         these_all_stamped_lick_times=[];
                                        %                                         these_all_stamped_lick_times=all_stamped_lick_times(these_per_corr,:);
                                        %                                         handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).stamped_lick_times=these_all_stamped_lick_times;
                                        %
                                        %                                         these_all_stamped_lick_ii=[];
                                        %                                         these_all_stamped_lick_ii=all_stamped_lick_ii(1,these_per_corr);
                                        %                                         handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).stamped_lick_ii=these_all_stamped_lick_ii;
                                        
                                        %Linear discriminant analysis for
                                        %trough PAC power
                                        
                                        these_all_log_P_timecoursePACwave=zeros(length(handles.drgbchoices.which_electrodes),sum(these_per_corr),length(t));
                                        these_all_log_P_timecoursePACwave(:,:,:)=all_log_P_timecoursePACwavetrough(PACii,:,these_per_corr,:);
                                        
                                        these_all_which_events=zeros(length(handles.drgbchoices.events_to_discriminate),sum(these_per_corr));
                                        for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                            kk=find(handles.drgbchoices.evTypeNos==handles.drgbchoices.events_to_discriminate(ii));
                                            these_all_which_events(ii,:)= all_which_eventsPACwave(kk,these_per_corr);
                                        end
                                        
                                        
                                        
                                        fprintf(1, ['LDA processed for trough PAC wavelet power for mouse No %d ' handles.drgbchoices.group_no_names{groupNo} ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' with %d trials\n'],mouseNo,N);
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
                                            
                                            %LFP power per trial per electrode
                                            measurements=zeros(N,length(handles.drgbchoices.which_electrodes));
                                            measurements(:,:)=these_all_log_P_timecoursePACwave(:,:,time_point)';
                                            
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
                                        
                                        if PACii==3
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
                                            
                                            %title(['LDA % correct for ' handles.drgbchoices.bwlabels{PACii} ' mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                                            xlabel('Time (sec)')
                                            ylabel('% correct trough')
                                            
                                            subplot(2,5,7)
                                            hold on
                                            
                                            plot(t,auROC,'-b')
                                            
                                            %Odor on markers
                                            plot([0 0],[0 0.5],'-k')
                                            odorhl=plot([0 2.5],[0 0],'-k','LineWidth',5);
                                            plot([2.5 2.5],[0 0.5],'-k')
                                            
                                            %title(['auROC for LDA for ' handles.drgbchoices.bwlabels{PACii} ' mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
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
                                            
                                            %title(['auROC for LDA for ' handles.drgbchoices.bwlabels{PACii} ' mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
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
                                            
                                            suptitle(['PAC wavelet power LDA analysis for Theta/' handles.drgbchoices.PACnames{PACii} ' PAC, mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                                        end
                                        
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).p_val_trough=p_val_trough;
                                        
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated=1;
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_trough=zeros(1,length(t));
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_trough(1,:)=discriminant_correct(1,:);
                                        
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).dimensionality_trough=zeros(1,length(t));
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).dimensionality_trough(1,:)=dimensionality(1,:);
                                        
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).auROC_trough=zeros(1,length(t));
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).auROC_trough=auROC;
                                        
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_shuffled_trough=zeros(1,length(t));
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_shuffled_trough(1,:)=discriminant_correct_shuffled(1,:);
                                        
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).test_out_per_timepoint_trough=test_out_per_timepoint;
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).shuffled_out_per_timepoint_trough=shuffled_out_per_timepoint;
                                        
                                        
                                        %                                         handles_out.t_power=t';
                                        %                                         handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).no_trials=N;
                                        %                                         handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).per_targets=per_targets;
                                        %
                                        %                                         these_all_which_events=[];
                                        %                                         these_all_which_events=all_which_eventsPACwave(:,these_per_corr);
                                        %                                         handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).which_events=these_all_which_events;
                                        %
                                        %                                         these_all_stamped_lick_times=[];
                                        %                                         these_all_stamped_lick_times=all_stamped_lick_times(these_per_corr,:);
                                        %                                         handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).stamped_lick_times=these_all_stamped_lick_times;
                                        %
                                        %                                         these_all_stamped_lick_ii=[];
                                        %                                         these_all_stamped_lick_ii=all_stamped_lick_ii(1,these_per_corr);
                                        %                                         handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).stamped_lick_ii=these_all_stamped_lick_ii;
                                        
                                        
                                    end
                                    
                                    if sum(handles.drgbchoices.which_discriminant==11)>0
                                        
                                        %PCA for peak PAC wavelet power
                                        
                                        these_all_log_P_timecoursePACwave=zeros(length(handles.drgbchoices.which_electrodes),sum(these_per_corr),length(t));
                                        these_all_log_P_timecoursePACwave(:,:,:)=all_log_P_timecoursePACwavepeak(PACii,:,these_per_corr,:);
                                        
                                        these_all_which_events=zeros(length(handles.drgbchoices.events_to_discriminate),sum(these_per_corr));
                                        for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                            kk=find(handles.drgbchoices.evTypeNos==handles.drgbchoices.events_to_discriminate(ii));
                                            these_all_which_events(ii,:)= all_which_eventsPACwave(kk,these_per_corr);
                                        end
                                        
                                        fprintf(1, ['PCA processed for peak PAC wavelet power for mouse No %d ' handles.drgbchoices.group_no_names{groupNo} ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' with %d trials\n'],mouseNo,N);
                                        fprintf(1,'For these events: ')
                                        for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                            fprintf(1,[handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(ii)} ' '])
                                        end
                                        fprintf(1,'\n')
                                        
                                        for time_point=1:length(t)
                                            par_t_out(time_point).principal_components=zeros(N,length(handles.drgbchoices.which_electrodes));
                                        end
                                        
                                        
                                        for time_point=1:length(t)
                                            
                                            %LFP power per trial per electrode
                                            measurements=zeros(N,length(handles.drgbchoices.which_electrodes));
                                            measurements(:,:)=these_all_log_P_timecoursePACwave(:,:,time_point)';
                                            
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
                                        
                                        principal_components=zeros(length(t),N,length(handles.drgbchoices.which_electrodes));
                                        PC_variance=zeros(length(t),length(handles.drgbchoices.which_electrodes));
                                        for time_point=1:length(t)
                                            principal_components(time_point,:,:)=par_t_out(time_point).principal_components;
                                            PC_variance(time_point,:)=par_t_out(time_point).PC_variance;
                                        end
                                        
                                        %Show the result of the PCA
                                        if PACii==3
                                            figNo=figNo+1;
                                            try
                                                close(figNo)
                                            catch
                                            end
                                            
                                            figure(figNo)
                                            
                                            %Show PCA before odor on
                                            subplot(2,4,5)
                                            hold on
                                            
                                            these_pcs=zeros(N,length(handles.drgbchoices.which_electrodes));
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
                                            
                                            these_pcs=zeros(N,length(handles.drgbchoices.which_electrodes));
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
                                        end
                                        
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PCA_calculated=1;
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).principal_components_peak=zeros(length(t),N,length(handles.drgbchoices.which_electrodes));
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).principal_components_peak=principal_components;
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).PC_variance_peak=zeros(length(t),length(handles.drgbchoices.which_electrodes));
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).PC_variance_peak=PC_variance;
                                        
                                        
                                        handles_out.t_power=t_power';
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).no_trials=N;
                                        
                                        these_all_which_events=[];
                                        these_all_which_events=all_which_eventsPACwave(:,these_per_corr);
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).which_events=these_all_which_events;
                                        
                                        %                                         these_all_stamped_lick_times=[];
                                        %                                         these_all_stamped_lick_times=all_stamped_lick_times(these_per_corr,:);
                                        %                                         handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).stamped_lick_times=these_all_stamped_lick_times;
                                        %
                                        %                                         these_all_stamped_lick_ii=[];
                                        %                                         these_all_stamped_lick_ii=all_stamped_lick_ii(1,these_per_corr);
                                        %                                         handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).stamped_lick_ii=these_all_stamped_lick_ii;
                                        
                                        
                                        %PCA for trough PAC wavelet power
                                        
                                        these_all_log_P_timecoursePACwave=zeros(length(handles.drgbchoices.which_electrodes),sum(these_per_corr),length(t));
                                        these_all_log_P_timecoursePACwave(:,:,:)=all_log_P_timecoursePACwavetrough(PACii,:,these_per_corr,:);
                                        
                                        these_all_which_events=zeros(length(handles.drgbchoices.events_to_discriminate),sum(these_per_corr));
                                        for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                            kk=find(handles.drgbchoices.evTypeNos==handles.drgbchoices.events_to_discriminate(ii));
                                            these_all_which_events(ii,:)= all_which_eventsPACwave(kk,these_per_corr);
                                        end
                                        
                                        fprintf(1, ['PCA processed for trough PAC wavelet power for mouse No %d ' handles.drgbchoices.group_no_names{groupNo} ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' with %d trials\n'],mouseNo,N);
                                        fprintf(1,'For these events: ')
                                        for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                            fprintf(1,[handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(ii)} ' '])
                                        end
                                        fprintf(1,'\n')
                                        
                                        for time_point=1:length(t)
                                            par_t_out(time_point).principal_components=zeros(N,length(handles.drgbchoices.which_electrodes));
                                        end
                                        
                                        for time_point=1:length(t)
                                            par_t_out(time_point).principal_components=[];
                                        end
                                        
                                        for time_point=1:length(t)
                                            
                                            %LFP power per trial per electrode
                                            measurements=zeros(N,length(handles.drgbchoices.which_electrodes));
                                            measurements(:,:)=these_all_log_P_timecoursePACwave(:,:,time_point)';
                                            
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
                                        
                                        principal_components=zeros(length(t),N,length(handles.drgbchoices.which_electrodes));
                                        PC_variance=zeros(length(t),length(handles.drgbchoices.which_electrodes));
                                        for time_point=1:length(t)
                                            principal_components(time_point,:,:)=par_t_out(time_point).principal_components;
                                            PC_variance(time_point,:)=par_t_out(time_point).PC_variance;
                                        end
                                        
                                        
                                        if PACii==3
                                            %Show PCA before odor on
                                            subplot(2,4,7)
                                            hold on
                                            
                                            these_pcs=zeros(N,length(handles.drgbchoices.which_electrodes));
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
                                            
                                            these_pcs=zeros(N,length(handles.drgbchoices.which_electrodes));
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
                                            suptitle(['PAC wavelet power PC1 for Theta/' handles.drgbchoices.PACnames{PACii} ' PAC, for mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                                            title('PC1 for trough')
                                            xlabel('Time (sec)')
                                            ylabel('PC1')
                                        end
                                        
                                        %                                         handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated=1;
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).principal_components_trough=zeros(length(t),N,length(handles.drgbchoices.which_electrodes));
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).principal_components_trough=principal_components;
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).PC_variance_trough=zeros(length(t),length(handles.drgbchoices.which_electrodes));
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).PC_variance_trough=PC_variance;
                                        
                                        %                                         handles_out.t_power=t';
                                        %                                         handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).no_trials=N;
                                        %
                                        %                                         these_all_which_events=[];
                                        %                                         these_all_which_events=all_which_eventsPACwave(:,these_per_corr);
                                        %                                         handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).which_events=these_all_which_events;
                                        %
                                        %                                         these_all_stamped_lick_times=[];
                                        %                                         these_all_stamped_lick_times=all_stamped_lick_times(these_per_corr,:);
                                        %                                         handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).stamped_lick_times=these_all_stamped_lick_times;
                                        %
                                        %                                         these_all_stamped_lick_ii=[];
                                        %                                         these_all_stamped_lick_ii=all_stamped_lick_ii(1,these_per_corr);
                                        %                                         handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).stamped_lick_ii=these_all_stamped_lick_ii;
                                        
                                    end
                                    
                                    if sum(handles.drgbchoices.which_discriminant==12)>0
                                        
                                        %First do peak PAC wavelet power
                                        %LDA for subsets of electrodes for PAC wavelet power
                                        
                                        t_from=2.15;
                                        t_to=2.5;
                                        subt=t((t>=t_from)&(t<=t_to));
                                        
                                        these_all_log_P_timecoursePACwave=zeros(length(handles.drgbchoices.which_electrodes),sum(these_per_corr),length(subt));
                                        these_all_log_P_timecoursePACwave(:,:,:)=all_log_P_timecoursePACwavepeak(PACii,:,these_per_corr,(t>=t_from)&(t<=t_to));
                                        
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
                                            
                                            suptitle(['PAC peak wavelet power LDA ' handles.drgbchoices.PACnames{PACii} ' mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                                        end
                                        
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).ecomb_discriminant_calculated_peak=1;
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).ecomb_par_out_peak=par_out;
                                        
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
                                            
                                            suptitle(['PAC wavelet power LDA ' handles.drgbchoices.PACnames{PACii} ' mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                                        end
                                        
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).ecomb_discriminant_calculated_trough=1;
                                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).ecomb_par_out_trough=par_out;
                                    end
                                    
                                else
                                    handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated=0;
                                    fprintf(1, ['LDA/PCA not processed for PAC wavelet power for mouse No %d ' handles.drgbchoices.group_no_names{groupNo} ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' because there were only %d trials (fewer than 20 trials)\n'],mouseNo,N);
                                end
                            end
                        end
                    end
                    
                    %Save all mouse data
                    if (sum(handles.drgbchoices.which_discriminant==15)>0)||(sum(handles.drgbchoices.which_discriminant==16)>0)
                        all_mouse_wav(mouseNo).group(groupNo).mouse_has_data=mouse_has_data;
                        if mouse_has_data==1
                            all_mouse_wav(mouseNo).group(groupNo).t=t_pac;
                            all_mouse_wav(mouseNo).group(groupNo).all_perCorr_pertrPACwave=all_perCorr_pertrPACwave;
                            all_mouse_wav(mouseNo).group(groupNo).all_licks_per_tPACwave=all_licks_per_tPACwave;
                            all_mouse_wav(mouseNo).group(groupNo).all_stamped_lick_iiPACwave=all_stamped_lick_iiPACwave;
                            all_mouse_wav(mouseNo).group(groupNo).all_log_P_timecoursePACwavepeak=all_log_P_timecoursePACwavepeak;
                            all_mouse_wav(mouseNo).group(groupNo).all_log_P_timecoursePACwavetrough=all_log_P_timecoursePACwavetrough;
                            all_mouse_wav(mouseNo).group(groupNo).all_which_eventsPACwave=all_which_eventsPACwave;
                        end
                    end
                    
                    %Calculate discriminant analysis and PCA for PAC
                    if (sum(handles.drgbchoices.which_discriminant==4)>0)||(sum(handles.drgbchoices.which_discriminant==5)>0)
                        
                        for PACii=1:length(handles.drgbchoices.PACburstLowF)
                            
                            for percent_correct_ii=1:length(handles.drgbchoices.per_lab)
                                
                                discriminant_correctPAC=zeros(1,length(t_pac));
                                discriminant_correct_shuffledPAC=zeros(1,length(t_pac));
                                auROCPAC=zeros(1,length(t_pac));
                                
                                these_per_corr=(all_perCorr_pertrPAC>=handles.drgbchoices.percent_windows(percent_correct_ii,1))...
                                    &(all_perCorr_pertrPAC<=handles.drgbchoices.percent_windows(percent_correct_ii,2));
                                N=sum(these_per_corr);
                                
                                %Do the analysis only if there are more than 20 trials
                                if N>=20
                                    
                                    
                                    if sum(handles.drgbchoices.which_discriminant==4)>0
                                        %Linear discriminant analysis for
                                        %phase
                                        
                                        these_all_phase_histo_timecourse=zeros(length(handles.drgbchoices.which_electrodes)*no_bins,sum(these_per_corr),length(t_pac));
                                        these_all_phase_histo_timecourse(:,:,:)=all_phase_histo_timecourse(PACii,:,these_per_corr,:);
                                        
                                        these_all_which_events=zeros(length(handles.drgbchoices.events_to_discriminate),sum(these_per_corr));
                                        for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                            kk=find(handles.drgbchoices.evTypeNos==handles.drgbchoices.events_to_discriminate(ii));
                                            these_all_which_events(ii,:)= all_which_eventsPAC(kk,these_per_corr);
                                        end
                                        
                                        test_out_per_timepoint=zeros(length(handles.drgbchoices.events_to_discriminate),N,length(t_pac));
                                        shuffled_out_per_timepoint=zeros(length(handles.drgbchoices.events_to_discriminate),N,length(t_pac));
                                        
                                        fprintf(1, ['LDA processed for phase for mouse No %d ' handles.drgbchoices.group_no_names{groupNo} ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' with %d trials\n'],mouseNo,N);
                                        fprintf(1,'For these events: ')
                                        for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                            fprintf(1,[handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(ii)} ' '])
                                        end
                                        fprintf(1,'\n')
                                        
                                        
                                        
                                        auROCPAC=zeros(1,length(t_pac));
                                        discriminant_correctPAC=zeros(1,length(t_pac));
                                        discriminant_correct_shuffledPAC=zeros(1,length(t_pac));
                                        test_out_per_timepoint=zeros(length(handles.drgbchoices.events_to_discriminate),N,length(t_pac));
                                        shuffled_out_per_timepoint=zeros(length(handles.drgbchoices.events_to_discriminate),N,length(t_pac));
                                        per_targets=zeros(length(handles.drgbchoices.events_to_discriminate),N);
                                        
                                        parfor time_point=1:length(t_pac)
                                            
                                            %LFP power per trial per electrode
                                            measurements=zeros(N,length(handles.drgbchoices.which_electrodes)*no_bins);
                                            measurements(:,:)=these_all_phase_histo_timecourse(:,:,time_point)';
                                            
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
                                                shuffled_measurements=zeros(N,length(handles.drgbchoices.which_electrodes)*no_bins);
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
                                            auROCPAC(1,time_point)=AUC-0.5;
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
                                            discriminant_correctPAC(1,time_point)=100*sum(sum(test_out.*per_targets))/N;
                                            discriminant_correct_shuffledPAC(1,time_point)=100*sum(sum(shuffled_out.*per_targets))/N;
                                            fprintf(1, 'LDA percent correct classification %d (for timepoint %d out of %d)\n',100*sum(per_targets(1,:)==test_out(1,:))/N,time_point,length(t_pac));
                                            pffft=1;
                                        end
                                    end
                                    
                                    
                                    
                                    if (sum(handles.drgbchoices.which_discriminant==4)>0)
                                        
                                        if PACii==3
                                            figNo=figNo+1;
                                            try
                                                close(figNo)
                                            catch
                                            end
                                            
                                            figure(figNo)
                                            
                                            subplot(1,2,1)
                                            hold on
                                            
                                            per95=prctile(discriminant_correct_shuffledPAC(1,:),95);
                                            per5=prctile(discriminant_correct_shuffledPAC(1,:),5);
                                            CIsh=[mean(discriminant_correct_shuffledPAC(1,:))-per5 per95-mean(discriminant_correct_shuffledPAC(1,:))]';
                                            [hlCR, hpCR] = boundedline([t_pac(1) t_pac(end)],[mean(discriminant_correct_shuffledPAC(1,:)) mean(discriminant_correct_shuffledPAC(1,:))], CIsh', 'r');
                                            
                                            plot(t_pac',discriminant_correctPAC(1,:),'-k')
                                            
                                            %Odor on markers
                                            plot([0 0],[0 100],'-k')
                                            odorhl=plot([0 2.5],[20 20],'-k','LineWidth',5);
                                            plot([2.5 2.5],[0 100],'-k')
                                            
                                            %title(['LDA % correct for ' handles.drgbchoices.bwlabels{bwii} ' mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                                            xlabel('Time (sec)')
                                            ylabel('Percent correct')
                                            
                                            subplot(1,2,2)
                                            hold on
                                            
                                            plot(t_pac,auROCPAC)
                                            
                                            %Odor on markers
                                            plot([0 0],[-0.3 0.5],'-k')
                                            odorhl=plot([0 2.5],[0 0],'-k','LineWidth',5);
                                            plot([2.5 2.5],[-0.3 0.5],'-k')
                                            
                                            %title(['auROC for phase LDA for ' handles.drgbchoices.bwlabels{bwii} ' mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                                            xlabel('Time (sec)')
                                            ylabel('auROC')
                                            ylim([-0.3 0.6])
                                            
                                            suptitle(['LDA phase analysis for ' handles.drgbchoices.PACnames{PACii} ' PAC mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                                        end
                                        
                                        handles_out.discriminant_per_mousePAC(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated=1;
                                        handles_out.discriminant_per_mousePAC(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correctPAC=zeros(1,length(t_pac));
                                        handles_out.discriminant_per_mousePAC(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correctPAC(1,:)=discriminant_correctPAC(1,:);
                                        
                                        handles_out.discriminant_per_mousePAC(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).auROC=zeros(1,length(t_pac));
                                        handles_out.discriminant_per_mousePAC(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).auROC=auROCPAC;
                                        
                                        handles_out.discriminant_per_mousePAC(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_shuffledPAC=zeros(1,length(t_pac));
                                        handles_out.discriminant_per_mousePAC(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_shuffledPAC(1,:)=discriminant_correct_shuffledPAC(1,:);
                                        
                                        handles_out.discriminant_per_mousePAC(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).test_out_per_timepoint=test_out_per_timepoint;
                                        handles_out.discriminant_per_mousePAC(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).shuffled_out_per_timepoint=shuffled_out_per_timepoint;
                                        
                                        handles_out.t_pac=t_pac';
                                        handles_out.discriminant_per_mousePAC(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).no_trials=N;
                                        handles_out.discriminant_per_mousePAC(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).per_targets=per_targets;
                                        
                                        these_all_which_events=[];
                                        these_all_which_events=all_which_eventsPAC(:,these_per_corr);
                                        handles_out.discriminant_per_mousePAC(mouseNo).group(groupNo).percent_correct(percent_correct_ii).which_events=these_all_which_events;
                                        
                                        these_all_stamped_lick_times=[];
                                        these_all_stamped_lick_times=all_stamped_lick_timesPAC(these_per_corr,:);
                                        handles_out.discriminant_per_mousePAC(mouseNo).group(groupNo).percent_correct(percent_correct_ii).stamped_lick_times=these_all_stamped_lick_times;
                                        
                                        these_all_stamped_lick_ii=[];
                                        these_all_stamped_lick_ii=all_stamped_lick_iiPAC(1,these_per_corr);
                                        handles_out.discriminant_per_mousePAC(mouseNo).group(groupNo).percent_correct(percent_correct_ii).stamped_lick_ii=these_all_stamped_lick_ii;
                                    end
                                    
                                    if sum(handles.drgbchoices.which_discriminant==5)>0
                                        
                                        %PCA
                                        
                                        these_all_phase_histo_timecourse=zeros(length(handles.drgbchoices.which_electrodes)*no_bins,sum(these_per_corr),length(t_pac));
                                        these_all_phase_histo_timecourse(:,:,:)=all_phase_histo_timecourse(PACii,:,these_per_corr,:);
                                        
                                        these_all_which_events=zeros(length(handles.drgbchoices.events_to_discriminate),sum(these_per_corr));
                                        for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                            kk=find(handles.drgbchoices.evTypeNos==handles.drgbchoices.events_to_discriminate(ii));
                                            these_all_which_events(ii,:)= all_which_eventsPAC(kk,these_per_corr);
                                        end
                                        
                                        test_out_per_timepoint=zeros(length(handles.drgbchoices.events_to_discriminate),N,length(t_pac));
                                        shuffled_out_per_timepoint=zeros(length(handles.drgbchoices.events_to_discriminate),N,length(t_pac));
                                        
                                        fprintf(1, ['PCA processed for PAC for mouse No %d ' handles.drgbchoices.group_no_names{groupNo} ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' with %d trials\n'],mouseNo,N);
                                        fprintf(1,'For these events: ')
                                        for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                            fprintf(1,[handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(ii)} ' '])
                                        end
                                        fprintf(1,'\n')
                                        
                                        for time_point=1:length(t_pac)
                                            par_t_out(time_point).principal_componentsPAC=zeros(N,length(handles.drgbchoices.which_electrodes)*no_bins);
                                        end
                                        
                                        
                                        for time_point=1:length(t_pac)
                                            par_t_out(time_point).principal_componentsPAC=[];
                                        end
                                        
                                        parfor time_point=1:length(t_pac)
                                            
                                            %LFP power per trial per electrode
                                            measurements=zeros(N,length(handles.drgbchoices.which_electrodes)*no_bins);
                                            measurements(:,:)=these_all_phase_histo_timecourse(:,:,time_point)';
                                            
                                            %Do the PCA
                                            [coeff,par_t_out(time_point).principal_componentsPAC,latent]=pca(measurements);
                                            
                                        end
                                        
                                        
                                        %Show a figure of the PCA and record the output
                                        szpto=size(par_t_out(time_point).principal_componentsPAC);
                                        principal_componentsPAC=zeros(length(t_pac),N,szpto(2));
                                        for time_point=1:length(t_pac)
                                            principal_componentsPAC(time_point,:,:)=par_t_out(time_point).principal_componentsPAC;
                                        end
                                        
                                        if PACii==3
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
                                            
                                            these_pcs=zeros(N,szpto(2));
                                            these_pcs(:,:)=principal_componentsPAC(6,:,:);
                                            plot(these_pcs(logical(these_all_which_events(1,:)),1),these_pcs(logical(these_all_which_events(1,:)),2),'or')
                                            plot(these_pcs(logical(these_all_which_events(2,:)),1),these_pcs(logical(these_all_which_events(2,:)),2),'ob')
                                            legend(handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(1)},handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(2)})
                                            xlabel('PC1')
                                            ylabel('PC2')
                                            title('-1 sec')
                                            
                                            %Show PCA after odor
                                            subplot(2,2,4)
                                            hold on
                                            
                                            these_pcs=zeros(N,szpto(2));
                                            these_pcs(:,:)=principal_componentsPAC(41,:,:);
                                            plot(these_pcs(logical(these_all_which_events(1,:)),1),these_pcs(logical(these_all_which_events(1,:)),2),'or')
                                            plot(these_pcs(logical(these_all_which_events(2,:)),1),these_pcs(logical(these_all_which_events(2,:)),2),'ob')
                                            legend(handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(1)},handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(2)})
                                            xlabel('PC1')
                                            ylabel('PC2')
                                            title('2.5 sec')
                                            
                                            %Show the timecourse for PC1
                                            subplot(2,2,[1,2])
                                            hold on
                                            
                                            PC1ev1=zeros(length(t_pac),sum(these_all_which_events(1,:)));
                                            PC1ev1(:,:)=principal_componentsPAC(:,logical(these_all_which_events(1,:)),1);
                                            
                                            PC1ev2=zeros(length(t_pac),sum(these_all_which_events(2,:)));
                                            PC1ev2(:,:)=principal_componentsPAC(:,logical(these_all_which_events(2,:)),1);
                                            
                                            mean_PC1ev2=mean(PC1ev2,2)';
                                            CIPC1ev2 = bootci(1000, {@mean, PC1ev2'});
                                            maxCIPC1ev2=max(CIPC1ev2(:));
                                            minCIPC1ev2=min(CIPC1ev2(:));
                                            CIPC1ev2(1,:)=mean_PC1ev2-CIPC1ev2(1,:);
                                            CIPC1ev2(2,:)=CIPC1ev2(2,:)-mean_PC1ev2;
                                            [hlCR, hpCR] = boundedline(t_pac',mean_PC1ev2', CIPC1ev2', 'b');
                                            
                                            mean_PC1ev1=mean(PC1ev1,2)';
                                            CIPC1ev1 = bootci(1000, {@mean, PC1ev1'});
                                            maxCIPC1ev1=max(CIPC1ev1(:));
                                            minCIPC1ev1=min(CIPC1ev1(:));
                                            CIPC1ev1(1,:)=mean_PC1ev1-CIPC1ev1(1,:);
                                            CIPC1ev1(2,:)=CIPC1ev1(2,:)-mean_PC1ev1;
                                            [hlCR, hpCR] = boundedline(t_pac',mean_PC1ev1', CIPC1ev1', 'r');
                                            
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
                                            title(['PC1 for ' handles.drgbchoices.PACnames{PACii} ' PAC mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                                            xlabel('Time (sec)')
                                            ylabel('PC1')
                                        end
                                        
                                        handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated=1;
                                        handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).principal_componentsPAC=zeros(length(t_pac),N,szpto(2));
                                        handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).principal_componentsPAC=principal_componentsPAC;
                                        
                                        handles_out.t_pac=t_pac';
                                        handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).no_trials=N;
                                        
                                        these_all_which_events=[];
                                        these_all_which_events=all_which_eventsPAC(:,these_per_corr);
                                        handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).which_events=these_all_which_events;
                                        
                                        these_all_stamped_lick_times=[];
                                        these_all_stamped_lick_times=all_stamped_lick_timesPAC(these_per_corr,:);
                                        handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).stamped_lick_times=these_all_stamped_lick_times;
                                        
                                        these_all_stamped_lick_ii=[];
                                        these_all_stamped_lick_ii=all_stamped_lick_iiPAC(1,these_per_corr);
                                        handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).stamped_lick_ii=these_all_stamped_lick_ii;
                                        
                                    end
                                    
                                    
                                    
                                else
                                    handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated=0;
                                    fprintf(1, ['LDA/PCA not processed for mouse No %d ' handles.drgbchoices.group_no_names{groupNo} ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' because there were only %d trials (fewer than 20 trials)\n'],mouseNo,N);
                                end
                            end
                        end
                    end
                    
                else
                    handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated=0;
                end
                save([handles.drgb.outPathName 'Discriminant_' handles.drgb.outFileName],'handles_out','-v7.3')
            end
        end
        
        
        if (mouse_has_data==1)&(process_data_for_mouse==1)& ((sum(handles.drgbchoices.which_discriminant==13)>0)||(sum(handles.drgbchoices.which_discriminant==14)>0))
            
            N_included=22;
            
            if (sum(all_groups_pertrPACwave==1)>=N_included)&(sum(all_groups_pertrPACwave==2)>=N_included)
                %Calculate discriminant analysis for wavelet power calculated at different PAC phases
                
                %Firts let's find performance of the LDA when the last N_included of group
                %1 and the first N_included trials of group 2 are probed independently
                
                t=t_pac;
                which_analysis_ii=1;
                for PACii=1:length(handles.drgbchoices.PACburstLowF)
                    
                    %Note, group 1 no laser, group 2 laser
                    for groupNo=1:2
                        %Group No 1 is before laser, group 2 is with laser light on
                        
                        discriminant_correct=zeros(1,length(t));
                        discriminant_correct_shuffled=zeros(1,length(t));
                        auROC=zeros(1,length(t));
                        
                        
                        trials_for_this_group=zeros(1,length(all_groups_pertrPACwave));
                        if groupNo==1
                            ii_last=find(all_groups_pertrPACwave==1,1,'last');
                            trials_for_this_group(1,ii_last-N_included+1:ii_last)=1;
                        else
                            ii_first=find(all_groups_pertrPACwave==2,1,'first');
                            trials_for_this_group(1,ii_first:N_included+ii_first-1)=1;
                        end
                        trials_for_this_group=logical(trials_for_this_group);
                        N=sum(trials_for_this_group);
                        
                        
                        %                     %Stamp the licks
                        %                     these_all_licks_per_tPACwave=[];
                        %                     these_all_licks_per_tPACwave=all_licks_per_tPACwave(trials_for_this_group,:);
                        %
                        %                     these_all_stamped_lick_ii=[];
                        %                     these_all_stamped_lick_ii=all_stamped_lick_iiPACwave(1,trials_for_this_group);
                        
                        
                        %Do the analysis only if there are more than 20 trials
                        if N>=20
                            
                            if (sum(handles.drgbchoices.which_discriminant==13)>0)
                                %Linear discriminant analysis for
                                %peak PAC power
                                
                                these_all_log_P_timecoursePACwave=zeros(length(handles.drgbchoices.which_electrodes),sum(trials_for_this_group),length(t));
                                these_all_log_P_timecoursePACwave(:,:,:)=all_log_P_timecoursePACwavepeak(PACii,:,trials_for_this_group,:);
                                
                                these_all_which_events=zeros(length(handles.drgbchoices.events_to_discriminate),sum(trials_for_this_group));
                                for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                    kk=find(handles.drgbchoices.evTypeNos==handles.drgbchoices.events_to_discriminate(ii));
                                    these_all_which_events(ii,:)= all_which_eventsPACwave(kk,trials_for_this_group);
                                end
                                
                                
                                if groupNo==1
                                    fprintf(1, ['LDA processed before laser for peak PAC wavelet power for mouse No %d ' handles.drgbchoices.group_no_names{groupNo} ' with %d trials\n'],mouseNo,N);
                                else
                                    fprintf(1, ['LDA processed during laser on for peak PAC wavelet power for mouse No %d ' handles.drgbchoices.group_no_names{groupNo}  ' with %d trials\n'],mouseNo,N);
                                end
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
                                    
                                    %LFP power per trial per electrode
                                    measurements=zeros(N,length(handles.drgbchoices.which_electrodes));
                                    measurements(:,:)=these_all_log_P_timecoursePACwave(:,:,time_point)';
                                    
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
                                    fprintf(1, 'LDA for peak PAC wavelet power percent correct classification %d (for timepoint %d out of %d)\n',100*sum(per_targets(1,:)==test_out(1,:))/N,time_point,length(t));
                                    
                                end
                                
                                
                                if PACii==3
                                    figNo=figNo+1;
                                    try
                                        close(figNo)
                                    catch
                                    end
                                    
                                    hFig=figure(figNo)
                                    set(hFig, 'units','normalized','position',[.1 .4 .6 .47])
                                    
                                    subplot(2,3,1)
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
                                    
                                    subplot(2,3,2)
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
                                    
                                    
                                    subplot(2,3,3)
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
                                end
                                
                                %Calculate p value for the wavelet power
                                p_val_peak=zeros(1,length(t));
                                for ii_t=1:length(t)
                                    splus_out=zeros(1,sum(these_all_which_events(1,:)==1));
                                    splus_out(1,:)=test_out_per_timepoint(1,these_all_which_events(1,:)==1,ii_t);
                                    sminus_out=zeros(1,sum(these_all_which_events(1,:)==0));
                                    sminus_out(1,:)=test_out_per_timepoint(1,these_all_which_events(1,:)==0,ii_t);
                                    p_val_peak(ii_t)=ranksum(splus_out,sminus_out);
                                end
                                
                                
                                %                             %Calculate p value for the licks
                                %                             p_val_lick=zeros(1,length(t));
                                %                             for ii_t=1:length(t)
                                %                                 splus_licks=zeros(1,sum((these_all_which_events(1,:)==1)&(these_all_stamped_lick_ii>0)));
                                %                                 splus_licks(1,:)=these_all_licks_per_tPACwave((these_all_which_events(1,:)==1)&(these_all_stamped_lick_ii>0),ii_t);
                                %                                 sminus_licks=zeros(1,sum((these_all_which_events(1,:)==0)&(these_all_stamped_lick_ii>0)));
                                %                                 sminus_licks(1,:)=these_all_licks_per_tPACwave((these_all_which_events(1,:)==0)&(these_all_stamped_lick_ii>0),ii_t);
                                %                                 p_val_lick(ii_t)=ranksum(splus_licks,sminus_licks);
                                %                             end
                                
                                
                                %suptitle(['PAC power LDA analysis for Theta/' handles.drgbchoices.PACnames{PACii} ' PAC, mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                                %                             handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).all_stamped_lick_ii=these_all_stamped_lick_ii;
                                %                             handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).all_licks_per_tPACwave=these_all_licks_per_tPACwave;
                                handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).p_val_peak=p_val_peak;
                                %                             handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).p_val_lick=p_val_lick;
                                
                                handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).which_analysis(which_analysis_ii).discriminant_calculated=1;
                                handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).discriminant_correct_peak=zeros(1,length(t));
                                handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).discriminant_correct_peak(1,:)=discriminant_correct(1,:);
                                
                                handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).dimensionality_peak=zeros(1,length(t));
                                handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).dimensionality_peak(1,:)=dimensionality(1,:);
                                
                                
                                handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).auROC_peak=zeros(1,length(t));
                                handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).auROC_peak=auROC;
                                
                                handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).discriminant_correct_shuffled_peak=zeros(1,length(t));
                                handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).discriminant_correct_shuffled_peak(1,:)=discriminant_correct_shuffled(1,:);
                                
                                handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).test_out_per_timepoint_peak=test_out_per_timepoint;
                                handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).shuffled_out_per_timepoint_peak=shuffled_out_per_timepoint;
                                
                                handles_out.t_power=t';
                                handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).no_trials=N;
                                handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).per_targets=per_targets;
                                
                                these_all_which_events=[];
                                these_all_which_events=all_which_eventsPACwave(:,trials_for_this_group);
                                handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).which_events=these_all_which_events;
                                
                                %                                         these_all_stamped_lick_times=[];
                                %                                         these_all_stamped_lick_times=all_stamped_lick_times(trials_for_this_group,:);
                                %                                         handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).stamped_lick_times=these_all_stamped_lick_times;
                                %
                                %                                         these_all_stamped_lick_ii=[];
                                %                                         these_all_stamped_lick_ii=all_stamped_lick_ii(1,trials_for_this_group);
                                %                                         handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).stamped_lick_ii=these_all_stamped_lick_ii;
                                
                                %Linear discriminant analysis for
                                %trough PAC power
                                
                                these_all_log_P_timecoursePACwave=zeros(length(handles.drgbchoices.which_electrodes),sum(trials_for_this_group),length(t));
                                these_all_log_P_timecoursePACwave(:,:,:)=all_log_P_timecoursePACwavetrough(PACii,:,trials_for_this_group,:);
                                
                                these_all_which_events=zeros(length(handles.drgbchoices.events_to_discriminate),sum(trials_for_this_group));
                                for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                    kk=find(handles.drgbchoices.evTypeNos==handles.drgbchoices.events_to_discriminate(ii));
                                    these_all_which_events(ii,:)= all_which_eventsPACwave(kk,trials_for_this_group);
                                end
                                
                                
                                
                                fprintf(1, ['LDA processed for trough PAC wavelet power for mouse No %d ' handles.drgbchoices.group_no_names{groupNo} ' ' handles.drgbchoices.per_lab{which_analysis_ii} ' with %d trials\n'],mouseNo,N);
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
                                    
                                    %LFP power per trial per electrode
                                    measurements=zeros(N,length(handles.drgbchoices.which_electrodes));
                                    measurements(:,:)=these_all_log_P_timecoursePACwave(:,:,time_point)';
                                    
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
                                
                                if PACii==3
                                    subplot(2,3,4)
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
                                    
                                    %title(['LDA % correct for ' handles.drgbchoices.bwlabels{PACii} ' mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{which_analysis_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                                    xlabel('Time (sec)')
                                    ylabel('% correct trough')
                                    
                                    subplot(2,3,5)
                                    hold on
                                    
                                    plot(t,auROC,'-b')
                                    
                                    %Odor on markers
                                    plot([0 0],[0 0.5],'-k')
                                    odorhl=plot([0 2.5],[0 0],'-k','LineWidth',5);
                                    plot([2.5 2.5],[0 0.5],'-k')
                                    
                                    %title(['auROC for LDA for ' handles.drgbchoices.bwlabels{PACii} ' mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{which_analysis_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                                    xlabel('Time (sec)')
                                    ylabel('auROC trough')
                                    ylim([-0.3 0.6])
                                    
                                    %Plot dimensionality
                                    subplot(2,3,6)
                                    hold on
                                    
                                    plot(t,dimensionality,'-b')
                                    
                                    mindim=min(dimensionality);
                                    maxdim=max(dimensionality);
                                    
                                    %Odor on markers
                                    plot([0 0],[mindim-0.1*(maxdim-mindim) maxdim+0.1*(maxdim-mindim)],'-k')
                                    odorhl=plot([0 2.5],[0 0],'-k','LineWidth',5);
                                    plot([2.5 2.5],[mindim-0.1*(maxdim-mindim) maxdim+0.1*(maxdim-mindim)],'-k')
                                    
                                    %title(['auROC for LDA for ' handles.drgbchoices.bwlabels{PACii} ' mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{which_analysis_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                                    xlabel('Time (sec)')
                                    ylabel('Dimensionality trough')
                                    xlim([-2 5])
                                    
                                    %                                 %Plot p value
                                    %                                 subplot(2,5,[4,5,9,10])
                                    %                                 hold on
                                    %
                                    %                                 p1=plot(t,log10(p_val_lick),'-k');
                                    %                                 p2=plot(t,log10(p_val_trough),'-b');
                                    %                                 p3=plot(t,log10(p_val_peak),'-m');
                                    %                                 plot([t(1) t(end)],[log10(0.05) log10(0.05)],'-r')
                                    %                                 legend([p1 p2 p3],{'Licks','Trough','Peak'})
                                    %                                 ylabel('log(p)')
                                    %                                 xlabel('Time (sec)')
                                    %                                 xlim([-2 5])
                                    
                                    if groupNo==1
                                        suptitle(['PAC wavelet power LDA analysis before laser for Theta/' handles.drgbchoices.PACnames{PACii} ' PAC, mouse No ' num2str(mouseNo) ' '  handles.drgbchoices.group_no_names{groupNo}])
                                    else
                                        suptitle(['PAC wavelet power LDA analysis during laser for Theta/' handles.drgbchoices.PACnames{PACii} ' PAC, mouse No ' num2str(mouseNo) ' '  handles.drgbchoices.group_no_names{groupNo}])
                                        
                                    end
                                end
                                
                                handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).p_val_trough=p_val_trough;
                                
                                handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).discriminant_calculated=1;
                                handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).discriminant_correct_trough=zeros(1,length(t));
                                handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).discriminant_correct_trough(1,:)=discriminant_correct(1,:);
                                
                                handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).dimensionality_trough=zeros(1,length(t));
                                handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).dimensionality_trough(1,:)=dimensionality(1,:);
                                
                                handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).auROC_trough=zeros(1,length(t));
                                handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).auROC_trough=auROC;
                                
                                handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).discriminant_correct_shuffled_trough=zeros(1,length(t));
                                handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).discriminant_correct_shuffled_trough(1,:)=discriminant_correct_shuffled(1,:);
                                
                                handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).test_out_per_timepoint_trough=test_out_per_timepoint;
                                handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).shuffled_out_per_timepoint_trough=shuffled_out_per_timepoint;
                                
                                
                                %                                         handles_out.t_power=t';
                                %                                         handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).no_trials=N;
                                %                                         handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).per_targets=per_targets;
                                %
                                %                                         these_all_which_events=[];
                                %                                         these_all_which_events=all_which_eventsPACwave(:,trials_for_this_group);
                                %                                         handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).which_events=these_all_which_events;
                                %
                                %                                         these_all_stamped_lick_times=[];
                                %                                         these_all_stamped_lick_times=all_stamped_lick_times(trials_for_this_group,:);
                                %                                         handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).stamped_lick_times=these_all_stamped_lick_times;
                                %
                                %                                         these_all_stamped_lick_ii=[];
                                %                                         these_all_stamped_lick_ii=all_stamped_lick_ii(1,trials_for_this_group);
                                %                                         handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).stamped_lick_ii=these_all_stamped_lick_ii;
                                
                                
                            end
                            
                            if sum(handles.drgbchoices.which_discriminant==14)>0
                                
                                %PCA for peak PAC wavelet power
                                
                                these_all_log_P_timecoursePACwave=zeros(length(handles.drgbchoices.which_electrodes),sum(trials_for_this_group),length(t));
                                these_all_log_P_timecoursePACwave(:,:,:)=all_log_P_timecoursePACwavepeak(PACii,:,trials_for_this_group,:);
                                
                                these_all_which_events=zeros(length(handles.drgbchoices.events_to_discriminate),sum(trials_for_this_group));
                                for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                    kk=find(handles.drgbchoices.evTypeNos==handles.drgbchoices.events_to_discriminate(ii));
                                    these_all_which_events(ii,:)= all_which_eventsPACwave(kk,trials_for_this_group);
                                end
                                
                                fprintf(1, ['PCA processed for peak PAC wavelet power for mouse No %d ' handles.drgbchoices.group_no_names{groupNo} ' ' handles.drgbchoices.per_lab{which_analysis_ii} ' with %d trials\n'],mouseNo,N);
                                fprintf(1,'For these events: ')
                                for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                    fprintf(1,[handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(ii)} ' '])
                                end
                                fprintf(1,'\n')
                                
                                for time_point=1:length(t)
                                    par_t_out(time_point).principal_components=zeros(N,length(handles.drgbchoices.which_electrodes));
                                end
                                
                                
                                for time_point=1:length(t)
                                    
                                    %LFP power per trial per electrode
                                    measurements=zeros(N,length(handles.drgbchoices.which_electrodes));
                                    measurements(:,:)=these_all_log_P_timecoursePACwave(:,:,time_point)';
                                    
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
                                
                                principal_components=zeros(length(t),N,length(handles.drgbchoices.which_electrodes));
                                PC_variance=zeros(length(t),length(handles.drgbchoices.which_electrodes));
                                for time_point=1:length(t)
                                    principal_components(time_point,:,:)=par_t_out(time_point).principal_components;
                                    PC_variance(time_point,:)=par_t_out(time_point).PC_variance;
                                end
                                
                                %Show the result of the PCA
                                if PACii==3
                                    figNo=figNo+1;
                                    try
                                        close(figNo)
                                    catch
                                    end
                                    
                                    figure(figNo)
                                    
                                    %Show PCA before odor on
                                    subplot(2,4,5)
                                    hold on
                                    
                                    these_pcs=zeros(N,length(handles.drgbchoices.which_electrodes));
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
                                    
                                    these_pcs=zeros(N,length(handles.drgbchoices.which_electrodes));
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
                                end
                                
                                handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PCA_calculated=1;
                                handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).principal_components_peak=zeros(length(t),N,length(handles.drgbchoices.which_electrodes));
                                handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).principal_components_peak=principal_components;
                                handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).PC_variance_peak=zeros(length(t),length(handles.drgbchoices.which_electrodes));
                                handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).PC_variance_peak=PC_variance;
                                
                                
                                handles_out.t_power=t_power';
                                handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).no_trials=N;
                                
                                these_all_which_events=[];
                                these_all_which_events=all_which_eventsPACwave(:,trials_for_this_group);
                                handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).which_events=these_all_which_events;
                                
                                %                                         these_all_stamped_lick_times=[];
                                %                                         these_all_stamped_lick_times=all_stamped_lick_times(trials_for_this_group,:);
                                %                                         handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).stamped_lick_times=these_all_stamped_lick_times;
                                %
                                %                                         these_all_stamped_lick_ii=[];
                                %                                         these_all_stamped_lick_ii=all_stamped_lick_ii(1,trials_for_this_group);
                                %                                         handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).stamped_lick_ii=these_all_stamped_lick_ii;
                                
                                
                                %PCA for trough PAC wavelet power
                                
                                these_all_log_P_timecoursePACwave=zeros(length(handles.drgbchoices.which_electrodes),sum(trials_for_this_group),length(t));
                                these_all_log_P_timecoursePACwave(:,:,:)=all_log_P_timecoursePACwavetrough(PACii,:,trials_for_this_group,:);
                                
                                these_all_which_events=zeros(length(handles.drgbchoices.events_to_discriminate),sum(trials_for_this_group));
                                for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                    kk=find(handles.drgbchoices.evTypeNos==handles.drgbchoices.events_to_discriminate(ii));
                                    these_all_which_events(ii,:)= all_which_eventsPACwave(kk,trials_for_this_group);
                                end
                                
                                fprintf(1, ['PCA processed for trough PAC wavelet power for mouse No %d ' handles.drgbchoices.group_no_names{groupNo} ' ' handles.drgbchoices.per_lab{which_analysis_ii} ' with %d trials\n'],mouseNo,N);
                                fprintf(1,'For these events: ')
                                for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                    fprintf(1,[handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(ii)} ' '])
                                end
                                fprintf(1,'\n')
                                
                                for time_point=1:length(t)
                                    par_t_out(time_point).principal_components=zeros(N,length(handles.drgbchoices.which_electrodes));
                                end
                                
                                for time_point=1:length(t)
                                    par_t_out(time_point).principal_components=[];
                                end
                                
                                for time_point=1:length(t)
                                    
                                    %LFP power per trial per electrode
                                    measurements=zeros(N,length(handles.drgbchoices.which_electrodes));
                                    measurements(:,:)=these_all_log_P_timecoursePACwave(:,:,time_point)';
                                    
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
                                
                                principal_components=zeros(length(t),N,length(handles.drgbchoices.which_electrodes));
                                PC_variance=zeros(length(t),length(handles.drgbchoices.which_electrodes));
                                for time_point=1:length(t)
                                    principal_components(time_point,:,:)=par_t_out(time_point).principal_components;
                                    PC_variance(time_point,:)=par_t_out(time_point).PC_variance;
                                end
                                
                                
                                if PACii==3
                                    %Show PCA before odor on
                                    subplot(2,4,7)
                                    hold on
                                    
                                    these_pcs=zeros(N,length(handles.drgbchoices.which_electrodes));
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
                                    
                                    these_pcs=zeros(N,length(handles.drgbchoices.which_electrodes));
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
                                    if groupNo==1
                                        suptitle(['PAC wavelet power before laser PC1 for Theta/' handles.drgbchoices.PACnames{PACii} ' PAC, for mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.group_no_names{groupNo}])
                                    else
                                        suptitle(['PAC wavelet power during laser PC1 for Theta/' handles.drgbchoices.PACnames{PACii} ' PAC, for mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.group_no_names{groupNo}])
                                    end
                                    title('PC1 for trough')
                                    xlabel('Time (sec)')
                                    ylabel('PC1')
                                end
                                
                                %                                         handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).discriminant_calculated=1;
                                handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).principal_components_trough=zeros(length(t),N,length(handles.drgbchoices.which_electrodes));
                                handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).principal_components_trough=principal_components;
                                handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).PC_variance_trough=zeros(length(t),length(handles.drgbchoices.which_electrodes));
                                handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).PC_variance_trough=PC_variance;
                                
                                %                                         handles_out.t_power=t';
                                %                                         handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).no_trials=N;
                                %
                                %                                         these_all_which_events=[];
                                %                                         these_all_which_events=all_which_eventsPACwave(:,trials_for_this_group);
                                %                                         handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).which_events=these_all_which_events;
                                %
                                %                                         these_all_stamped_lick_times=[];
                                %                                         these_all_stamped_lick_times=all_stamped_lick_times(trials_for_this_group,:);
                                %                                         handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).stamped_lick_times=these_all_stamped_lick_times;
                                %
                                %                                         these_all_stamped_lick_ii=[];
                                %                                         these_all_stamped_lick_ii=all_stamped_lick_ii(1,trials_for_this_group);
                                %                                         handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).stamped_lick_ii=these_all_stamped_lick_ii;
                                
                            end
                            
                            
                            
                        else
                            handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).discriminant_calculated=0;
                            fprintf(1, ['LDA/PCA not processed for PAC wavelet power for mouse No %d ' handles.drgbchoices.group_no_names{groupNo} ' ' handles.drgbchoices.per_lab{which_analysis_ii} ' because there were only %d trials (fewer than 20 trials)\n'],mouseNo,N);
                        end
                    end
                end
                
                %Now do the analysis where you train the before laser and ask what the
                %prediction is for the after laser
                which_analysis_ii=2;
                for PACii=1:length(handles.drgbchoices.PACburstLowF)
                    
                    %Train on group 1, predict groupNo 2
                    groupNo=1;
                    
                    %Group No 1 is before laser, group 2 is with laser light on
                    
                    discriminant_correct=zeros(1,length(t));
                    discriminant_correct_shuffled=zeros(1,length(t));
                    auROC=zeros(1,length(t));
                    
                    trials_for_this_group=zeros(1,length(all_groups_pertrPACwave));
                    
                    ii_last=find(all_groups_pertrPACwave==1,1,'last');
                    trials_for_this_group(1,ii_last-N_included+1:ii_last)=1;
                    trials_for_this_group=logical(trials_for_this_group);
                    
                    N=sum(trials_for_this_group);
                    
                    trials_for_this_test_group=zeros(1,length(all_groups_pertrPACwave));
                    
                    ii_first=find(all_groups_pertrPACwave==2,1,'first');
                    trials_for_this_test_group(1,ii_first:N_included+ii_first-1)=1;
                    trials_for_this_test_group=logical(trials_for_this_test_group);
                    
                    
                    %Do the analysis only if there are more than 20 trials
                    if N>=20
                        
                        if (sum(handles.drgbchoices.which_discriminant==13)>0)
                            %Linear discriminant analysis for
                            %peak PAC power
                            
                            these_all_log_P_timecoursePACwave=zeros(length(handles.drgbchoices.which_electrodes),sum(trials_for_this_group),length(t));
                            these_all_log_P_timecoursePACwave(:,:,:)=all_log_P_timecoursePACwavepeak(PACii,:,trials_for_this_group,:);
                            
                            these_all_which_events=zeros(length(handles.drgbchoices.events_to_discriminate),sum(trials_for_this_group));
                            for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                kk=find(handles.drgbchoices.evTypeNos==handles.drgbchoices.events_to_discriminate(ii));
                                these_all_which_events(ii,:)= all_which_eventsPACwave(kk,trials_for_this_group);
                            end
                            
                            
                            these_all_test_log_P_timecoursePACwave=zeros(length(handles.drgbchoices.which_electrodes),sum(trials_for_this_test_group),length(t));
                            these_all_test_log_P_timecoursePACwave(:,:,:)=all_log_P_timecoursePACwavepeak(PACii,:,trials_for_this_test_group,:);
                            
                            
                            these_all_test_which_events=zeros(length(handles.drgbchoices.events_to_discriminate),sum(trials_for_this_test_group));
                            for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                kk=find(handles.drgbchoices.evTypeNos==handles.drgbchoices.events_to_discriminate(ii));
                                these_all_test_which_events(ii,:)= all_which_eventsPACwave(kk,trials_for_this_test_group);
                            end
                            
                            
                            
                            fprintf(1, ['LDA trained on trials before laser, predicted during laser for peak PAC wavelet power for mouse No %d ' handles.drgbchoices.group_no_names{groupNo} ' with %d trials\n'],mouseNo,N);
                            
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
                                
                                %LFP power per trial per electrode
                                measurements=zeros(N,length(handles.drgbchoices.which_electrodes));
                                measurements(:,:)=these_all_log_P_timecoursePACwave(:,:,time_point)';
                                
                                test_measurements=zeros(N,length(handles.drgbchoices.which_electrodes));
                                test_measurements(:,:)=these_all_test_log_P_timecoursePACwave(:,:,time_point)';
                                
                                %                                 %Dimensionality
                                %                                 %Rows: trials, Columns: electrodes
                                %                                 Signal=measurements;
                                %                                 dimensionality(time_point) = nansum(eig(cov(Signal)))^2/nansum(eig(cov(Signal)).^2);
                                
                                %Enter strings labeling each event (one event for
                                %each trial)
                                events=[];
                                test_events=[];
                                for ii=1:N
                                    this_event=zeros(length(handles.drgbchoices.events_to_discriminate),1);
                                    this_event=these_all_which_events(:,ii);
                                    
                                    events{ii,1}=handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(logical(this_event))};
                                    
                                    this_event=zeros(length(handles.drgbchoices.events_to_discriminate),1);
                                    this_event=these_all_test_which_events(:,ii);
                                    
                                    test_events{ii,1}=handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(logical(this_event))};
                                    
                                end
                                
                                tested_events=[];
                                %                             shuffled_tested_events=[];
                                scores=[];
                                
                                %Train with all pre laser
                                idxTrn=ones(N,1);
                                
                                %Store the training data in a table.
                                tblTrn=[];
                                tblTrn = array2table(measurements(logical(idxTrn),:));
                                tblTrn.Y = events(logical(idxTrn));
                                
                                %Train a discriminant analysis model using the training set and default options.
                                %By default this is a regularized linear discriminant analysis (LDA)
                                Mdl = fitcdiscr(tblTrn,'Y');
                                
                                for ii=1:N
                                    
                                    %Test the laser events
                                    idxTest=zeros(N,1);
                                    idxTest(ii)=1;
                                    
                                    %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
                                    [label,score] = predict(Mdl,test_measurements(logical(idxTest),:));
                                    
                                    tested_events{ii,1}=label{1};
                                    scores(ii)=score(2);
                                    
                                end
                                
                                %Calculate auROC
                                [X,Y,T,AUC] = perfcurve(test_events,scores',handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(2)});
                                auROC(1,time_point)=AUC-0.5;
                                %test_out_per_timepoint=zeros(length(handles.drgbchoices.events_to_discriminate),N,length(t));
                                
                                
                                per_targets=these_all_test_which_events;
                                
                                test_out=zeros(length(handles.drgbchoices.events_to_discriminate),N);
                                %                                 shuffled_out=zeros(length(handles.drgbchoices.events_to_discriminate),N);
                                for ii=1:N
                                    for jj=1:length(handles.drgbchoices.events_to_discriminate)
                                        if strcmp(handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(jj)},tested_events{ii})
                                            test_out(jj,ii)=1;
                                        else
                                            test_out(jj,ii)=0;
                                        end
                                        %                                         if strcmp(handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(jj)},shuffled_tested_events{ii})
                                        %                                             shuffled_out(jj,ii)=1;
                                        %                                         else
                                        %                                             shuffled_out(jj,ii)=0;
                                        %                                         end
                                    end
                                end
                                
                                test_out_per_timepoint(:,:,time_point)=test_out;
                                %                                 shuffled_out_per_timepoint(:,:,time_point)=shuffled_out;
                                discriminant_correct(1,time_point)=100*sum(sum(test_out.*per_targets))/N;
                                %                                 discriminant_correct_shuffled(1,time_point)=100*sum(sum(shuffled_out.*per_targets))/N;
                                fprintf(1, 'LDA for peak PAC wavelet power percent correct classification %d (for timepoint %d out of %d)\n',100*sum(per_targets(1,:)==test_out(1,:))/N,time_point,length(t));
                                
                            end
                            
                            
                            if PACii==3
                                figNo=figNo+1;
                                try
                                    close(figNo)
                                catch
                                end
                                
                                hFig=figure(figNo)
                                set(hFig, 'units','normalized','position',[.1 .4 .6 .47])
                                
                                subplot(2,2,1)
                                hold on
                                
                                %                                 per95=prctile(discriminant_correct_shuffled(1,:),95);
                                %                                 per5=prctile(discriminant_correct_shuffled(1,:),5);
                                %                                 CIsh=[mean(discriminant_correct_shuffled(1,:))-per5 per95-mean(discriminant_correct_shuffled(1,:))]';
                                %                                 [hlCR, hpCR] = boundedline([t(1) t(end)],[mean(discriminant_correct_shuffled(1,:)) mean(discriminant_correct_shuffled(1,:))], CIsh', 'r');
                                %
                                plot(t',discriminant_correct(1,:),'-k')
                                
                                %Odor on markers
                                plot([0 0],[0 100],'-k')
                                odorhl=plot([0 2.5],[20 20],'-k','LineWidth',5);
                                plot([2.5 2.5],[0 100],'-k')
                                
                                %title(['LDA % correct for ' handles.drgbchoices.bwlabels{PACii} ' mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                                xlabel('Time (sec)')
                                ylabel('% correct peak')
                                
                                subplot(2,2,2)
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
                                
                                
                                %                             subplot(2,5,3)
                                %                             hold on
                                %
                                %                             plot(t,dimensionality,'-b')
                                %
                                %                             maxdim=max(dimensionality);
                                %                             mindim=min(dimensionality);
                                %
                                %                             %Odor on markers
                                %                             plot([0 0],[mindim-0.1*(maxdim-mindim) maxdim+0.1*(maxdim-mindim)],'-k')
                                %                             odorhl=plot([0 2.5],[mindim mindim],'-k','LineWidth',5);
                                %                             plot([2.5 2.5],[mindim-0.1*(maxdim-mindim) maxdim+0.1*(maxdim-mindim)],'-k')
                                %
                                %                             %title(['auROC for LDA for ' handles.drgbchoices.bwlabels{PACii} ' mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                                %                             xlabel('Time (sec)')
                                %                             ylabel('dimensionality peak')
                                %                             ylim([mindim-0.1*(maxdim-mindim) maxdim+0.1*(maxdim-mindim)])
                            end
                            
                            %Calculate p value for the wavelet power
                            p_val_peak=zeros(1,length(t));
                            for ii_t=1:length(t)
                                splus_out=zeros(1,sum(these_all_which_events(1,:)==1));
                                splus_out(1,:)=test_out_per_timepoint(1,these_all_which_events(1,:)==1,ii_t);
                                sminus_out=zeros(1,sum(these_all_which_events(1,:)==0));
                                sminus_out(1,:)=test_out_per_timepoint(1,these_all_which_events(1,:)==0,ii_t);
                                p_val_peak(ii_t)=ranksum(splus_out,sminus_out);
                            end
                            
                            
                            %                         %Calculate p value for the licks
                            %                         p_val_lick=zeros(1,length(t));
                            %                         for ii_t=1:length(t)
                            %                             splus_licks=zeros(1,sum((these_all_which_events(1,:)==1)&(these_all_stamped_lick_ii>0)));
                            %                             splus_licks(1,:)=these_all_licks_per_tPACwave((these_all_which_events(1,:)==1)&(these_all_stamped_lick_ii>0),ii_t);
                            %                             sminus_licks=zeros(1,sum((these_all_which_events(1,:)==0)&(these_all_stamped_lick_ii>0)));
                            %                             sminus_licks(1,:)=these_all_licks_per_tPACwave((these_all_which_events(1,:)==0)&(these_all_stamped_lick_ii>0),ii_t);
                            %                             p_val_lick(ii_t)=ranksum(splus_licks,sminus_licks);
                            %                         end
                            %
                            
                            %suptitle(['PAC power LDA analysis for Theta/' handles.drgbchoices.PACnames{PACii} ' PAC, mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                            %                         handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).all_stamped_lick_ii=these_all_stamped_lick_ii;
                            %                         handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).all_licks_per_tPACwave=these_all_licks_per_tPACwave;
                            handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).p_val_peak=p_val_peak;
                            %                         handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).p_val_lick=p_val_lick;
                            
                            handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).which_analysis(which_analysis_ii).discriminant_calculated=1;
                            handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).discriminant_correct_peak=zeros(1,length(t));
                            handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).discriminant_correct_peak(1,:)=discriminant_correct(1,:);
                            
                            %                             handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).dimensionality_peak=zeros(1,length(t));
                            %                             handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).dimensionality_peak(1,:)=dimensionality(1,:);
                            
                            
                            handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).auROC_peak=zeros(1,length(t));
                            handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).auROC_peak=auROC;
                            
                            %                             handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).discriminant_correct_shuffled_peak=zeros(1,length(t));
                            %                             handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).discriminant_correct_shuffled_peak(1,:)=discriminant_correct_shuffled(1,:);
                            %
                            handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).test_out_per_timepoint_peak=test_out_per_timepoint;
                            %                             handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).shuffled_out_per_timepoint_peak=shuffled_out_per_timepoint;
                            
                            handles_out.t_power=t';
                            handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).no_trials=N;
                            handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).per_targets=per_targets;
                            
                            these_all_which_events=[];
                            these_all_which_events=all_which_eventsPACwave(:,trials_for_this_group);
                            handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).which_events=these_all_test_which_events;
                            
                            %                                         these_all_stamped_lick_times=[];
                            %                                         these_all_stamped_lick_times=all_stamped_lick_times(trials_for_this_group,:);
                            %                                         handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).stamped_lick_times=these_all_stamped_lick_times;
                            %
                            %                                         these_all_stamped_lick_ii=[];
                            %                                         these_all_stamped_lick_ii=all_stamped_lick_ii(1,trials_for_this_group);
                            %                                         handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).stamped_lick_ii=these_all_stamped_lick_ii;
                            
                            %Linear discriminant analysis for
                            %Trough PAC power
                            
                            these_all_log_P_timecoursePACwave=zeros(length(handles.drgbchoices.which_electrodes),sum(trials_for_this_group),length(t));
                            these_all_log_P_timecoursePACwave(:,:,:)=all_log_P_timecoursePACwavetrough(PACii,:,trials_for_this_group,:);
                            
                            these_all_which_events=zeros(length(handles.drgbchoices.events_to_discriminate),sum(trials_for_this_group));
                            for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                kk=find(handles.drgbchoices.evTypeNos==handles.drgbchoices.events_to_discriminate(ii));
                                these_all_which_events(ii,:)= all_which_eventsPACwave(kk,trials_for_this_group);
                            end
                            
                            these_all_test_log_P_timecoursePACwave=zeros(length(handles.drgbchoices.which_electrodes),sum(trials_for_this_test_group),length(t));
                            these_all_test_log_P_timecoursePACwave(:,:,:)=all_log_P_timecoursePACwavetrough(PACii,:,trials_for_this_test_group,:);
                            
                            
                            these_all_test_which_events=zeros(length(handles.drgbchoices.events_to_discriminate),sum(trials_for_this_test_group));
                            for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                kk=find(handles.drgbchoices.evTypeNos==handles.drgbchoices.events_to_discriminate(ii));
                                these_all_test_which_events(ii,:)= all_which_eventsPACwave(kk,trials_for_this_test_group);
                            end
                            
                            
                            fprintf(1, ['LDA trained with trials before laser, prediction with trials during laser for trough PAC wavelet power for mouse No %d ' handles.drgbchoices.group_no_names{groupNo}  ' with %d trials\n'],mouseNo,N);
                            
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
                                
                                %LFP power per trial per electrode
                                measurements=zeros(N,length(handles.drgbchoices.which_electrodes));
                                measurements(:,:)=these_all_log_P_timecoursePACwave(:,:,time_point)';
                                
                                test_measurements=zeros(N,length(handles.drgbchoices.which_electrodes));
                                test_measurements(:,:)=these_all_test_log_P_timecoursePACwave(:,:,time_point)';
                                
                                %                                 %Dimensionality
                                %                                 %Rows: trials, Columns: electrodes
                                %                                 Signal=measurements;
                                %                                 dimensionality(time_point) = nansum(eig(cov(Signal)))^2/nansum(eig(cov(Signal)).^2);
                                
                                %Enter strings labeling each event (one event for
                                %each trial)
                                events=[];
                                test_events=[];
                                
                                for ii=1:N
                                    this_event=zeros(length(handles.drgbchoices.events_to_discriminate),1);
                                    this_event=these_all_which_events(:,ii);
                                    
                                    events{ii,1}=handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(logical(this_event))};
                                    
                                    this_event=zeros(length(handles.drgbchoices.events_to_discriminate),1);
                                    this_event=these_all_test_which_events(:,ii);
                                    
                                    test_events{ii,1}=handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(logical(this_event))};
                                    
                                end
                                
                                tested_events=[];
                                
                                scores=[];
                                
                                %Train with all pre laser
                                idxTrn=ones(N,1);
                                
                                %Store the training data in a table.
                                tblTrn=[];
                                tblTrn = array2table(measurements(logical(idxTrn),:));
                                tblTrn.Y = events(logical(idxTrn));
                                
                                %Train a discriminant analysis model using the training set and default options.
                                %By default this is a regularized linear discriminant analysis (LDA)
                                Mdl = fitcdiscr(tblTrn,'Y');
                                
                                for ii=1:N
                                    %Predict for trials during laser
                                    idxTest=zeros(N,1);
                                    idxTest(ii)=1;
                                    
                                    %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
                                    [label,score] = predict(Mdl,test_measurements(logical(idxTest),:));
                                    
                                    tested_events{ii,1}=label{1};
                                    scores(ii)=score(2);
                                    
                                end
                                
                                %Calculate auROC
                                [X,Y,T,AUC] = perfcurve(test_events,scores',handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(2)});
                                auROC(1,time_point)=AUC-0.5;
                                %test_out_per_timepoint=zeros(length(handles.drgbchoices.events_to_discriminate),N,length(t));
                                
                                
                                per_targets=these_all_test_which_events;
                                
                                test_out=zeros(length(handles.drgbchoices.events_to_discriminate),N);
                                
                                for ii=1:N
                                    for jj=1:length(handles.drgbchoices.events_to_discriminate)
                                        if strcmp(handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(jj)},tested_events{ii})
                                            test_out(jj,ii)=1;
                                        else
                                            test_out(jj,ii)=0;
                                        end
                                        
                                    end
                                end
                                
                                test_out_per_timepoint(:,:,time_point)=test_out;
                                
                                discriminant_correct(1,time_point)=100*sum(sum(test_out.*per_targets))/N;
                                
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
                            
                            if PACii==3
                                subplot(2,2,3)
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
                                
                                %title(['LDA % correct for ' handles.drgbchoices.bwlabels{PACii} ' mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{which_analysis_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                                xlabel('Time (sec)')
                                ylabel('% correct trough')
                                
                                subplot(2,2,4)
                                hold on
                                
                                plot(t,auROC,'-b')
                                
                                %Odor on markers
                                plot([0 0],[0 0.5],'-k')
                                odorhl=plot([0 2.5],[0 0],'-k','LineWidth',5);
                                plot([2.5 2.5],[0 0.5],'-k')
                                
                                %title(['auROC for LDA for ' handles.drgbchoices.bwlabels{PACii} ' mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{which_analysis_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                                xlabel('Time (sec)')
                                ylabel('auROC trough')
                                ylim([-0.3 0.6])
                                
                                
                                %
                                %                             %Plot p value
                                %                             subplot(2,5,[4,5,9,10])
                                %                             hold on
                                %
                                %                             p1=plot(t,log10(p_val_lick),'-k');
                                %                             p2=plot(t,log10(p_val_trough),'-b');
                                %                             p3=plot(t,log10(p_val_peak),'-m');
                                %                             plot([t(1) t(end)],[log10(0.05) log10(0.05)],'-r')
                                %                             legend([p1 p2 p3],{'Licks','Trough','Peak'})
                                %                             ylabel('log(p)')
                                %                             xlabel('Time (sec)')
                                %                             xlim([-2 5])
                                
                                
                                suptitle(['PAC wavelet power LDA analysis for during laser predictions (trained before laser) for Theta/' handles.drgbchoices.PACnames{PACii} ' PAC, mouse No ' num2str(mouseNo) ' '  handles.drgbchoices.group_no_names{groupNo}])
                                
                            end
                            
                            handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).p_val_trough=p_val_trough;
                            
                            handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).discriminant_calculated=1;
                            handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).discriminant_correct_trough=zeros(1,length(t));
                            handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).discriminant_correct_trough(1,:)=discriminant_correct(1,:);
                            
                            %                             handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).dimensionality_trough=zeros(1,length(t));
                            %                             handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).dimensionality_trough(1,:)=dimensionality(1,:);
                            
                            handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).auROC_trough=zeros(1,length(t));
                            handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).auROC_trough=auROC;
                            
                            handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).discriminant_correct_shuffled_trough=zeros(1,length(t));
                            handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).discriminant_correct_shuffled_trough(1,:)=discriminant_correct_shuffled(1,:);
                            
                            handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).test_out_per_timepoint_trough=test_out_per_timepoint;
                            %                             handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).shuffled_out_per_timepoint_trough=shuffled_out_per_timepoint;
                            
                            
                            %                                         handles_out.t_power=t';
                            %                                         handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).no_trials=N;
                            %                                         handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).PACii(PACii).per_targets=per_targets;
                            %
                            %                                         these_all_which_events=[];
                            %                                         these_all_which_events=all_which_eventsPACwave(:,trials_for_this_group);
                            %                                         handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).which_events=these_all_which_events;
                            %
                            %                                         these_all_stamped_lick_times=[];
                            %                                         these_all_stamped_lick_times=all_stamped_lick_times(trials_for_this_group,:);
                            %                                         handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).stamped_lick_times=these_all_stamped_lick_times;
                            %
                            %                                         these_all_stamped_lick_ii=[];
                            %                                         these_all_stamped_lick_ii=all_stamped_lick_ii(1,trials_for_this_group);
                            %                                         handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).stamped_lick_ii=these_all_stamped_lick_ii;
                            
                            
                        end
                        
                        
                        
                        
                        
                    else
                        handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).which_analysis(which_analysis_ii).discriminant_calculated=0;
                        fprintf(1, ['LDA/PCA not processed for PAC wavelet power for mouse No %d ' handles.drgbchoices.group_no_names{groupNo} ' ' handles.drgbchoices.per_lab{which_analysis_ii} ' because there were only %d trials (fewer than 20 trials)\n'],mouseNo,N);
                    end
                    
                end
            else
                fprintf(1, ['LDA was not processed for mouse No %d because the number of trials is smaller than %d trials\n'],mouseNo,N_included);
            end
        end
    end
    
    %Calculate discriminant analysis for wavelet power calculated at
    %different PAC phases for data for all mice
    if (sum(handles.drgbchoices.which_discriminant==15)>0)||(sum(handles.drgbchoices.which_discriminant==16)>0)
        t=t_pac;
        for groupNo=1:max(handles.drgbchoices.group_no)
            for PACii=1:length(handles.drgbchoices.PACburstLowF)
                
                for percent_correct_ii=1:length(handles.drgbchoices.per_lab)
                    
                    %First determine how many trials per mouse
                    N=0;
                    mouse_included=[];
                    max_mouse=[];
                    these_all_licks_per_tPACwave=[];
                    these_all_stamped_lick_ii=[];
                    these_all_which_events_licks=[];
                    N_licks=0;
                    for mouseNo=1:max(handles.drgbchoices.mouse_no)
                        if all_mouse_wav(mouseNo).group(groupNo).mouse_has_data==1
                            these_per_corr=(all_mouse_wav(mouseNo).group(groupNo).all_perCorr_pertrPACwave>=handles.drgbchoices.percent_windows(percent_correct_ii,1))...
                                &(all_mouse_wav(mouseNo).group(groupNo).all_perCorr_pertrPACwave<=handles.drgbchoices.percent_windows(percent_correct_ii,2));
                            if N==sum(these_per_corr)
                                max_mouse=mouseNo;
                            else
                                if sum(these_per_corr)>=20
                                    mouse_included(mouseNo)=1;
                                    if N<sum(these_per_corr)
                                        max_mouse=mouseNo;
                                        N=sum(these_per_corr);
                                    end
                                else
                                    mouse_included(mouseNo)=0;
                                end
                            end
                            
                            N_licks=N_licks+sum(these_per_corr);
                            
                        end
                    end
                    
                    
                    these_all_licks_per_tPACwave=zeros(N_licks,length(t));
                    these_all_stamped_lick_ii=zeros(1,N_licks);
                    these_all_which_events_licks=zeros(1,N_licks);
                    ii_licks=1;
                    
                    %Now find these_all_log_P_timecoursePACwave and licks
                    these_all_log_P_timecoursePACwave_peak=zeros(length(handles.drgbchoices.which_electrodes)*sum(mouse_included),N,length(t));
                    these_all_log_P_timecoursePACwave_trough=zeros(length(handles.drgbchoices.which_electrodes)*sum(mouse_included),N,length(t));
                    
                    
                    no_mice_included=0;
                    for mouseNo=1:max(handles.drgbchoices.mouse_no)
                        if max_mouse==mouseNo
                            %This is the mouse with the largest number of trials
                            no_mice_included=no_mice_included+1;
                            first_electrode=(no_mice_included-1)*length(handles.drgbchoices.which_electrodes)+1;
                            
                            these_per_corr=(all_mouse_wav(mouseNo).group(groupNo).all_perCorr_pertrPACwave>=handles.drgbchoices.percent_windows(percent_correct_ii,1))...
                                &(all_mouse_wav(mouseNo).group(groupNo).all_perCorr_pertrPACwave<=handles.drgbchoices.percent_windows(percent_correct_ii,2));
                            
                            these_all_log_P_timecoursePACwave_peak(first_electrode:first_electrode+length(handles.drgbchoices.which_electrodes)-1,:,:)=...
                                all_mouse_wav(mouseNo).group(groupNo).all_log_P_timecoursePACwavepeak(PACii,:,these_per_corr,:);
                            
                            these_all_log_P_timecoursePACwave_trough(first_electrode:first_electrode+length(handles.drgbchoices.which_electrodes)-1,:,:)=...
                                all_mouse_wav(mouseNo).group(groupNo).all_log_P_timecoursePACwavetrough(PACii,:,these_per_corr,:);
                            
                            these_all_which_events=zeros(length(handles.drgbchoices.events_to_discriminate),sum(these_per_corr));
                            for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                kk=find(handles.drgbchoices.evTypeNos==handles.drgbchoices.events_to_discriminate(ii));
                                these_all_which_events(ii,:)= all_mouse_wav(mouseNo).group(groupNo).all_which_eventsPACwave(kk,these_per_corr);
                            end
                            
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
                            
                            these_all_licks_per_tPACwave(ii_licks:ii_licks+ size(all_mouse_wav(mouseNo).group(groupNo).all_licks_per_tPACwave(these_per_corr,:),1)-1,:)=...
                                all_mouse_wav(mouseNo).group(groupNo).all_licks_per_tPACwave(these_per_corr,:);
                            
                            these_all_stamped_lick_ii(1,ii_licks:ii_licks+ length(all_mouse_wav(mouseNo).group(groupNo).all_stamped_lick_iiPACwave(1,these_per_corr))-1)=...
                                all_mouse_wav(mouseNo).group(groupNo).all_stamped_lick_iiPACwave(1,these_per_corr);
                            
                            these_all_which_events_licks(1,ii_licks:ii_licks+ length(all_mouse_wav(mouseNo).group(groupNo).all_stamped_lick_iiPACwave(1,these_per_corr))-1)=...
                                these_all_which_events(1,:);
                            
                            ii_licks=ii_licks+length(all_mouse_wav(mouseNo).group(groupNo).all_stamped_lick_iiPACwave(1,these_per_corr));
                        else
                            if mouse_included(mouseNo)==1
                                %Enter the existing trials
                                no_mice_included=no_mice_included+1;
                                first_electrode=(no_mice_included-1)*length(handles.drgbchoices.which_electrodes)+1;
                                
                                those_per_corr=(all_mouse_wav(mouseNo).group(groupNo).all_perCorr_pertrPACwave>=handles.drgbchoices.percent_windows(percent_correct_ii,1))...
                                    &(all_mouse_wav(mouseNo).group(groupNo).all_perCorr_pertrPACwave<=handles.drgbchoices.percent_windows(percent_correct_ii,2));
                                
                                those_all_log_P_timecoursePACwave_peak=zeros(length(handles.drgbchoices.which_electrodes),sum(those_per_corr),length(t));
                                those_all_log_P_timecoursePACwave_peak(:,:,:)=...
                                    all_mouse_wav(mouseNo).group(groupNo).all_log_P_timecoursePACwavepeak(PACii,:,those_per_corr,:);
                                
                                those_all_log_P_timecoursePACwave_trough=zeros(length(handles.drgbchoices.which_electrodes),sum(those_per_corr),length(t));
                                those_all_log_P_timecoursePACwave_trough(:,:,:)=...
                                    all_mouse_wav(mouseNo).group(groupNo).all_log_P_timecoursePACwavetrough(PACii,:,those_per_corr,:);
                                
                                those_all_which_events=zeros(length(handles.drgbchoices.events_to_discriminate),sum(those_per_corr));
                                
                                for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                    kk=find(handles.drgbchoices.evTypeNos==handles.drgbchoices.events_to_discriminate(ii));
                                    those_all_which_events(ii,:)= all_mouse_wav(mouseNo).group(groupNo).all_which_eventsPACwave(kk,those_per_corr);
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
                                
                                these_all_licks_per_tPACwave(ii_licks:ii_licks+ size(all_mouse_wav(mouseNo).group(groupNo).all_licks_per_tPACwave(those_per_corr,:),1)-1,:)=...
                                    all_mouse_wav(mouseNo).group(groupNo).all_licks_per_tPACwave(those_per_corr,:);
                                
                                these_all_stamped_lick_ii(1,ii_licks:ii_licks+ length(all_mouse_wav(mouseNo).group(groupNo).all_stamped_lick_iiPACwave(1,those_per_corr))-1)=...
                                    all_mouse_wav(mouseNo).group(groupNo).all_stamped_lick_iiPACwave(1,those_per_corr);
                                
                                these_all_which_events_licks(1,ii_licks:ii_licks+ length(all_mouse_wav(mouseNo).group(groupNo).all_stamped_lick_iiPACwave(1,those_per_corr))-1)=...
                                    those_all_which_events(1,:);
                                
                                ii_licks=ii_licks+length(all_mouse_wav(mouseNo).group(groupNo).all_stamped_lick_iiPACwave(1,those_per_corr));
                            end
                        end
                    end
                    
                    No_electrodes=first_electrode+length(handles.drgbchoices.which_electrodes)-1;
                    
                    
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
                            
%                             parfor time_point=1:length(t)
                            for time_point=1:length(t)
                                %
                                %LFP power per trial per electrode
                                measurements=zeros(N,No_electrodes);
                                measurements(:,:)=these_all_log_P_timecoursePACwave_peak(:,:,time_point)';
                                
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
                            
                            these_all_which_events=[];
                            these_all_which_events=all_which_eventsPACwave(:,these_per_corr);
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).which_events=these_all_which_events;
                            
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
                                
                                %LFP power per trial per electrode
                                measurements=zeros(N,length(handles.drgbchoices.which_electrodes));
                                measurements(:,:)=these_all_log_P_timecoursePACwave_trough(:,:,time_point)';
                                
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
                                par_t_out(time_point).principal_components=zeros(N,length(handles.drgbchoices.which_electrodes));
                            end
                            
                            
                            for time_point=1:length(t)
                                
                                %LFP power per trial per electrode
                                measurements=zeros(N,length(handles.drgbchoices.which_electrodes));
                                measurements(:,:)=these_all_log_P_timecoursePACwavepeak(:,:,time_point)';
                                
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
                            
                            principal_components=zeros(length(t),N,length(handles.drgbchoices.which_electrodes));
                            PC_variance=zeros(length(t),length(handles.drgbchoices.which_electrodes));
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
                            
                            these_pcs=zeros(N,length(handles.drgbchoices.which_electrodes));
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
                            
                            these_pcs=zeros(N,length(handles.drgbchoices.which_electrodes));
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
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).principal_components_peak=zeros(length(t),N,length(handles.drgbchoices.which_electrodes));
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).principal_components_peak=principal_components;
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).PC_variance_peak=zeros(length(t),length(handles.drgbchoices.which_electrodes));
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).PC_variance_peak=PC_variance;
                            
                            
                            handles_out.t_power=t_power';
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).no_trials=N;
                            
                            these_all_which_events=[];
                            these_all_which_events=all_which_eventsPACwave(:,these_per_corr);
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).which_events=these_all_which_events;
                            
                            %                                         these_all_stamped_lick_times=[];
                            %                                         these_all_stamped_lick_times=all_stamped_lick_times(these_per_corr,:);
                            %                                         handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).stamped_lick_times=these_all_stamped_lick_times;
                            %
                            %                                         these_all_stamped_lick_ii=[];
                            %                                         these_all_stamped_lick_ii=all_stamped_lick_ii(1,these_per_corr);
                            %                                         handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).stamped_lick_ii=these_all_stamped_lick_ii;
                            
                            
                            %PCA for trough PAC wavelet power
                            
                            
                            
                            fprintf(1, ['PCA processed for trough PAC wavelet power for all mice for ' handles.drgbchoices.group_no_names{groupNo} ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' with %d trials\n'],N);
                            fprintf(1,'For these events: ')
                            for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                fprintf(1,[handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(ii)} ' '])
                            end
                            fprintf(1,'\n')
                            
                            for time_point=1:length(t)
                                par_t_out(time_point).principal_components=zeros(N,length(handles.drgbchoices.which_electrodes));
                            end
                            
                            for time_point=1:length(t)
                                par_t_out(time_point).principal_components=[];
                            end
                            
                            for time_point=1:length(t)
                                
                                %LFP power per trial per electrode
                                measurements=zeros(N,length(handles.drgbchoices.which_electrodes));
                                measurements(:,:)=these_all_log_P_timecoursePACwavetrough(:,:,time_point)';
                                
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
                            
                            principal_components=zeros(length(t),N,length(handles.drgbchoices.which_electrodes));
                            PC_variance=zeros(length(t),length(handles.drgbchoices.which_electrodes));
                            for time_point=1:length(t)
                                principal_components(time_point,:,:)=par_t_out(time_point).principal_components;
                                PC_variance(time_point,:)=par_t_out(time_point).PC_variance;
                            end
                            
                            
                            
                            %Show PCA before odor on
                            subplot(2,4,7)
                            hold on
                            
                            these_pcs=zeros(N,length(handles.drgbchoices.which_electrodes));
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
                            
                            these_pcs=zeros(N,length(handles.drgbchoices.which_electrodes));
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
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).principal_components_trough=zeros(length(t),N,length(handles.drgbchoices.which_electrodes));
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).principal_components_trough=principal_components;
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).PC_variance_trough=zeros(length(t),length(handles.drgbchoices.which_electrodes));
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
                        
                    else
                        handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated=0;
                        fprintf(1, ['LDA/PCA not processed for PAC wavelet power for all mice ' handles.drgbchoices.group_no_names{groupNo} ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' because there were only %d trials (fewer than 20 trials)\n'],N);
                    end
                end
            end
        end
    end
    
end

pffft1=1;


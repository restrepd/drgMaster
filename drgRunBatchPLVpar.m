function drgRunBatchPLVpar



tic

first_file=1;

[choiceFileName,choiceBatchPathName] = uigetfile({'drgbChoicesPLV*.m'},'Select the .m file with all the choices for analysis');
fprintf(1, ['\ndrgRunBatchPLV run for ' choiceFileName '\n\n']);

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
    %     this_jt=handles.drgbchoices.FileName{filNum};
    %     handles.temp_exist(filNum)=exist([handles.drgb.outPathName tempDirName '/temp_' this_jt(10:end)]);
    %
    %     if handles.temp_exist(filNum)==2
    %         %If it was processed load the temp result
    %         load([handles.drgb.outPathName tempDirName '/temp_' this_jt(10:end)])
    %         lfp_per_file(filNum)=this_lfp_per_file;
    %
    %     end
end

%Now process each mouse separately
if all_files_present==1
    
    %Generates a PLV analysis for all electrode pairs
    odorOn=2;
    decimation_factor=40;
    handles.decimation_factor=decimation_factor;
    
    gcp
    
    
    parfor mouse_ii=1:max(handles.drgbchoices.mouse_no)
        %     for mouse_ii=1:max(handles.drgbchoices.mouse_no)
        if exist([handles.drgb.outPathName tempDirName '/' 'temp_mouse' num2str(mouse_ii) '.mat'],'file')~=2
            handlespf=struct();
            handlespf=handles;
            
            no_files=handlespf.drgbchoices.no_files;
            all_files=1:no_files;
            
            
            PLV_per_mouse=[];
            ii_PLV=0;
            
            handlespf.time_start=handlespf.drgbchoices.time_start;
            handlespf.time_end=handlespf.drgbchoices.time_end;
            
            
            for bw_ii=1:length(handlespf.drgbchoices.lowF)
                lowF1=handlespf.drgbchoices.lowF(bw_ii);
                lowF2=handlespf.drgbchoices.highF(bw_ii);
                for per_ii=1:size(handlespf.drgbchoices.percent_windows,1)
                    for eventNo_ii=1:length(handlespf.drgbchoices.evTypeNos)
                        
                        ii_vet=0;
                        these_LFPs=[];
                        found_trials=0;
                        
                        for filNum=first_file:no_files
                            
                            
                            if handlespf.drgbchoices.mouse_no(filNum)==mouse_ii
                                
                                
                                this_jt=handlespf.drgbchoices.FileName{filNum};
                                
                                
                                %Othrwise read the jt_times and do processing
                                %read the jt_times file
                                jtFileName=handlespf.drgbchoices.FileName{filNum};
                                if iscell(handlespf.drgbchoices.PathName)
                                    jtPathName=handlespf.drgbchoices.PathName{filNum};
                                else
                                    jtPathName=handlespf.drgbchoices.PathName;
                                end
                                
                                drgRead_jt_times(jtPathName,jtFileName);
                                FileName=[jtFileName(10:end-4) '_drg.mat'];
                                fullName=[jtPathName,FileName];
                                my_drg={'drg'};
                                S=load(fullName,my_drg{:});
                                handlespf.drg=S.drg;
                                
                                sessionNo=handlespf.sessionNo;
                                Fs=handlespf.drg.session(sessionNo).draq_p.ActualRate;
                                
                                %Enter trials
                                firstTr=1;
                                lastTr=handlespf.drg.drta_p.trialNo;
                                
                                [perCorr, encoding_trials, retrieval_trials, encoding_this_evTypeNo,retrieval_this_evTypeNo]=drgFindEncRetr(handlespf);
                                
                                
                                
                                
                                handlespf.evTypeNo=handlespf.drgbchoices.evTypeNos(eventNo_ii);
                                vetted_trNos=[];
                                vetted_evNos=[];
                                
                                
                                
                                for trNo=firstTr:lastTr
                                    pCorr=perCorr(drgFindEvNo(handlespf,trNo,sessionNo,odorOn));
                                    if (pCorr>=handlespf.drgbchoices.percent_windows(per_ii,1))&(pCorr<=handlespf.drgbchoices.percent_windows(per_ii,2))
                                        evNo = drgFindEvNo(handlespf,trNo,sessionNo);
                                        if evNo~=-1
                                            excludeTrial=drgExcludeTrialLFP(handlespf.drg,handlespf.peakLFPNo,handlespf.drg.session(sessionNo).events(handlespf.evTypeNo).times(evNo),sessionNo);
                                            
                                            if excludeTrial==0
                                                
                                                ii_elect=0;
                                                all_can_read=1;
                                                while (all_can_read==1)&(ii_elect<handlespf.drgbchoices.no_electrodes)
                                                    handlespf.peakLFPNo=ii_elect+1;
                                                    [LFP, trialNo, can_read] = drgGetTrialLFPData(handlespf, handlespf.peakLFPNo, evNo, handlespf.evTypeNo, handlespf.time_start, handlespf.time_end);
                                                    if (can_read==1)
                                                        ii_elect=ii_elect+1;
                                                        LFP=decimate(LFP,decimation_factor);
                                                        these_LFPs(ii_elect,1:length(LFP),ii_vet+1)=LFP;
                                                    else
                                                        all_can_read=0;
                                                    end
                                                end
                                                
                                                if all_can_read==1
                                                    if found_trials==0
                                                        ii_PLV=ii_PLV+1;
                                                        found_trials=1;
                                                        PLV_per_mouse.PLV(ii_PLV).eventNo_ii=eventNo_ii;
                                                        PLV_per_mouse.PLV(ii_PLV).per_ii=per_ii;
                                                        PLV_per_mouse.PLV(ii_PLV).bw_ii=bw_ii;
                                                        PLV_per_mouse.group_no=handles.drgbchoices.group_no(filNum);
                                                        PLV_per_mouse.mouse_no=mouse_ii;
                                                    end
                                                    
                                                    ii_vet=ii_vet+1;
                                                    PLV_per_mouse.PLV(ii_PLV).no_trials=ii_vet;
                                                    PLV_per_mouse.PLV(ii_PLV).perCorr(ii_vet)=perCorr(drgFindEvNo(handlespf,trNo,sessionNo,odorOn));
                                                    PLV_per_mouse.PLV(ii_PLV).trial(ii_vet).trialNo=trNo;
                                                    PLV_per_mouse.PLV(ii_PLV).trial(ii_vet).fileNo=filNum;
                                                    
                                                    vetted_trNos(ii_vet)=trNo;
                                                    vetted_evNos(ii_vet)=evNo;
                                                    if handlespf.displayData==1
                                                        this_file_per_mouse=find(all_files(handles.drgbchoices.mouse_no==mouse_ii)==filNum);
                                                        fprintf(1, ['Mouse %d, file number %d//%d, bandwidth %d, per_ii %d, evNo %d, trial %d\n'], mouse_ii, this_file_per_mouse, sum(handles.drgbchoices.mouse_no==mouse_ii), bw_ii, per_ii, eventNo_ii,ii_vet);
                                                    end
                                                end
                                                
                                            end
                                            
                                        end
                                    end
                                    %end %if eventstamps...
                                end %for trNo
                                
                                
                            end
                        end
                        
                        %Now add shuffled trials to LFP handlespf.drgbchoices.no_electrodes+1
                        if found_trials==1
                            perm_ii_vet = randperm(length(vetted_trNos));
                            ii_elec=1;
                            for jj_vet=1:length(vetted_trNos)
                                these_LFPs(handlespf.drgbchoices.no_electrodes+1,1:length(LFP),jj_vet)=these_LFPs(ii_elec,1:length(LFP),perm_ii_vet(jj_vet));
                                ii_elec=ii_elec+1;
                                if ii_elec>handlespf.drgbchoices.no_electrodes
                                    ii_elec=1;
                                end
                            end
                            
                            %                             if ~isempty(these_LFPs)
                            %Calculate PLV
                            filtSpec=struct();
                            filtSpec.range = [lowF1 lowF2];
                            filtSpec.order=50;
                            [plv,delta_phase]=pn_eegPLV(these_LFPs, Fs/decimation_factor, filtSpec);
                            
                            %Discard the ends
                            ii_start=floor(handlespf.time_pad*Fs/decimation_factor)+1;
                            ii_end=ii_start+floor((handlespf.time_end-handlespf.time_start-2*handlespf.time_pad)*Fs/decimation_factor)-1;
                            plv = plv(ii_start:ii_end,:,:);
                            delta_phase=delta_phase(ii_start:ii_end,:,:);
                            PLV_per_mouse.PLV(ii_PLV).plv=plv;
                            PLV_per_mouse.PLV(ii_PLV).delta_phase=delta_phase;
                            
                            time=handlespf.time_start+handlespf.time_pad:1/floor(Fs/decimation_factor):handlespf.time_end-handlespf.time_pad;
                            time=time(1:length(plv));
                            PLV_per_mouse.PLV(ii_PLV).time=time;
                        end
                    end
                    
                end
                
            end
            
            if handlespf.displayData==1
                fprintf(1, ['Finished processing mouse %d\n'], mouse_ii);
            end
            
            %Save this temp file
            drgSavePLVPar([handlespf.drgb.outPathName tempDirName '/'],['mouse' num2str(mouse_ii)],PLV_per_mouse,mouse_ii)
            PLV_per_mouse=[];
        end
    end
    
    handles_out=[];
    ii_PLV=0;
    for mouse_ii=1:max(handles.drgbchoices.mouse_no)
        load([handles.drgb.outPathName tempDirName '/' 'temp_mouse' num2str(mouse_ii) '.mat'])
        for jj=1:length(this_PLV_per_mouse.PLV)
            this_PLV_per_mouse.PLV(jj).mouse_no=mouse_ii;
            this_PLV_per_mouse.PLV(jj).group_no=this_PLV_per_mouse.group_no;
            ii_PLV=ii_PLV+1;
            handles_out.PLV(ii_PLV)=this_PLV_per_mouse.PLV(jj);
        end
    end
    
    handles_out.drg=handles.drg;
    save([handles.drgb.outPathName handles.drgb.outFileName],'handles_out','-v7.3')
    

    
end









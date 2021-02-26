function drgRunBatchPLV



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
    figNo=0;
    %Empty vectors
    handles_out.PLV=[];
    ii_PLV=0;
    
    for mouse_ii=1:max(handles.drgbchoices.mouse_no)
        
        
        no_files=handles.drgbchoices.no_files;
        
        for filNum=first_file:no_files
            
            
            if handles.drgbchoices.mouse_no(filNum)==mouse_ii

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
                
                sessionNo=handles.sessionNo;
                Fs=handles.drg.session(sessionNo).draq_p.ActualRate;
                
                %Enter trials
                firstTr=1;
                lastTr=handles.drg.drta_p.trialNo;
                
                [perCorr, encoding_trials, retrieval_trials, encoding_this_evTypeNo,retrieval_this_evTypeNo]=drgFindEncRetr(handles);
                
                this_evTypeNo=handles.evTypeNo;
                t_start=handles.time_start;
                t_end=handles.time_end;
                handles.time_start=handles.drgbchoices.time_start;
                handles.time_end=handles.drgbchoices.time_end;
                
                
                for bw_ii=1:length(handles.drgbchoices.lowF)
                    lowF1=handles.drgbchoices.lowF(bw_ii);
                    lowF2=handles.drgbchoices.highF(bw_ii);
                    for per_ii=1:size(handles.drgbchoices.percent_windows,1)
                        for eventNo_ii=1:length(handles.drgbchoices.evTypeNos)
                            ii_PLV=ii_PLV+1;
                            handles_out.PLV(ii_PLV).eventNo_ii=eventNo_ii;
                            handles_out.PLV(ii_PLV).per_ii=per_ii;
                            handles_out.PLV(ii_PLV).bw_ii=bw_ii;
                            handles.evTypeNo=handles.drgbchoices.evTypeNos(eventNo_ii);
                            vetted_trNos=[];
                            vetted_evNos=[];
                            ii_vet=0;
                            these_LFPs=[];
                            for trNo=firstTr:lastTr
                                pCorr=perCorr(drgFindEvNo(handles,trNo,sessionNo,odorOn));
                                if (pCorr>=handles.drgbchoices.percent_windows(per_ii,1))&(pCorr<=handles.drgbchoices.percent_windows(per_ii,2))
                                    evNo = drgFindEvNo(handles,trNo,sessionNo);
                                    if evNo~=-1
                                        excludeTrial=drgExcludeTrialLFP(handles.drg,handles.peakLFPNo,handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);
                                        
                                        if excludeTrial==0
                                            
                                            ii_elect=0;
                                            all_can_read=1;
                                            while (all_can_read==1)&(ii_elect<handles.drgbchoices.no_electrodes)
                                                handles.peakLFPNo=ii_elect+1;
                                                [LFP, trialNo, can_read] = drgGetTrialLFPData(handles, handles.peakLFPNo, evNo, handles.evTypeNo, handles.time_start, handles.time_end);
                                                if (can_read==1)
                                                    ii_elect=ii_elect+1;
                                                    LFP=decimate(LFP,decimation_factor);
                                                    these_LFPs(ii_elect,1:length(LFP),ii_vet+1)=LFP;
                                                else
                                                    all_can_read=0;
                                                end
                                            end
                                            
                                            if all_can_read==1
                                                ii_vet=ii_vet+1;
                                                handles_out.PLV(ii_PLV).no_trials=ii_vet;
                                                handles_out.PLV(ii_PLV).perCorr(ii_vet)=perCorr(drgFindEvNo(handles,trNo,sessionNo,odorOn));
                                                handles_out.PLV(ii_PLV).trial(ii_vet).trialNo=trNo;
                                                
                                                vetted_trNos(ii_vet)=trNo;
                                                vetted_evNos(ii_vet)=evNo;
                                                if handles.displayData==1
                                                    fprintf(1, ['Mouse %d, file %d, bandwidth %d, per_ii %d, evNo %d, trial %d\n'], mouse_ii, filNum, bw_ii, per_ii, eventNo_ii,ii_vet);
                                                end
                                            end
                                            
                                        end
                                        
                                    end
                                end
                                %end %if eventstamps...
                            end %for trNo
                            
                            %Now add shuffled trials to LFP handles.drgbchoices.no_electrodes+1
                            perm_ii_vet = randperm(length(vetted_trNos));
                            ii_elec=1;
                            for ii_vet=1:length(vetted_trNos)
                                these_LFPs(handles.drgbchoices.no_electrodes+1,1:length(LFP),ii_vet)=these_LFPs(ii_elec,1:length(LFP),perm_ii_vet(ii_vet));
                                ii_elec=ii_elec+1;
                                if ii_elec>handles.drgbchoices.no_electrodes
                                    ii_elec=1;
                                end
                            end
                            
                            if ~isempty(these_LFPs)
                                %Calculate PLV
                                filtSpec.range = [lowF1 lowF2];
                                filtSpec.order=50;
                                [plv,delta_phase]=pn_eegPLV(these_LFPs, Fs/decimation_factor, filtSpec);
                                
                                %Discard the ends
                                ii_start=floor(handles.time_pad*Fs/decimation_factor)+1;
                                ii_end=ii_start+floor((handles.time_end-handles.time_start-2*handles.time_pad)*Fs/decimation_factor)-1;
                                plv = plv(ii_start:ii_end,:,:);
                                delta_phase=delta_phase(ii_start:ii_end,:,:);
                                handles_out.PLV(ii_PLV).plv=plv;
                                handles_out.PLV(ii_PLV).delta_phase=delta_phase;
                                
                                time=handles.time_start+handles.time_pad:1/floor(Fs/decimation_factor):handles.time_end-handles.time_pad;
                                time=time(1:length(plv));
                                handles_out.PLV(ii_PLV).time=time;
                            end
                            
                        end
                    end
                end
                
                handles.evTypeNo=this_evTypeNo;
                handles.time_start=t_start;
                handles.time_end=t_end;
                
                %Save this temp file
                save([handles.drgb.outPathName handles.drgb.outFileName],'handles_out','-v7.3')

            end
            
        end
        
        %Display the data
        
        if handles.displayData==1
            
            
            
            %Plot PLV timecourse between brain regions
            
            
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            
            
            hFig = figure(figNo);
            set(hFig, 'units','normalized','position',[.05 .05 .45 .7])
            
            this_sub=0;
            
            for bw_ii=1:length(handles.drgbchoices.lowF)
                
                for per_ii=size(handles.drgbchoices.percent_windows,1):-1:1
                    this_sub=this_sub+1;
                    subplot(length(handles.drgbchoices.lowF),size(handles.drgbchoices.percent_windows,1),this_sub)
                    hold on
                    for eventNo_ii=length(handles.drgbchoices.evTypeNos):-1:1
                        
                        %Get the PLVs
                        for ii=1:ii_PLV
                            if (handles_out.PLV(ii).eventNo_ii==eventNo_ii)&(handles_out.PLV(ii).per_ii==per_ii)&(handles_out.PLV(ii).bw_ii==bw_ii)
                                this_ii_PLV=ii;
                            end
                        end
                        
                        
                        these_PLVs=zeros(length(handles.drgbchoices.reference_electrodes)*length(handles.drgbchoices.other_electrodes),length(handles_out.PLV(this_ii_PLV).time));
                        ii_PLVs=0;
                        for ii_ref=handles.drgbchoices.reference_electrodes
                            for ii_oth=handles.drgbchoices.other_electrodes
                                this_plv=zeros(1,length(handles_out.PLV(this_ii_PLV).time));
                                if ii_ref<ii_oth
                                    this_plv(1,:)=handles_out.PLV(this_ii_PLV).plv(:,ii_ref,ii_oth);
                                else
                                    this_plv(1,:)=handles_out.PLV(this_ii_PLV).plv(:,ii_oth,ii_ref);
                                end
                                ii_PLVs=ii_PLVs+1;
                                these_PLVs(ii_PLVs,:)=this_plv;
                            end
                        end
                        
                        %Calculate plv for 100 msec intervals
                        dt=0.1;
                        no_time_points=(Fs/decimation_factor)*dt;
                        these_PLVs_dec=zeros(size(these_PLVs,1),size(these_PLVs,2)/no_time_points);
                        time_dec=zeros(1,size(these_PLVs,2)/no_time_points);
                        time=handles_out.PLV(this_ii_PLV).time;
                        
                        for elec_pair=1:size(these_PLVs,1)
                            for ii_dec_tp=1:size(these_PLVs,2)/no_time_points
                                these_PLVs_dec(elec_pair,ii_dec_tp)=mean(these_PLVs(elec_pair,1+(ii_dec_tp-1)*no_time_points:ii_dec_tp*no_time_points));
                                if elec_pair==1
                                    time_dec(1,ii_dec_tp)=mean(time(1+(ii_dec_tp-1)*no_time_points:ii_dec_tp*no_time_points));
                                end
                            end
                        end
                        
                        
                        mean_plv=mean(these_PLVs_dec);
                        
                        
                        if eventNo_ii==2
                            plot(time_dec,mean_plv, 'b');
                        else
                            plot(time_dec,mean_plv, 'r');
                        end
                        
                        
                        %Get PLS
                        these_PLVs=zeros(length(handles.drgbchoices.reference_electrodes),length(handles_out.PLV(this_ii_PLV).time));
                        ii_PLVs=0;
                        for ii_ref=1:handles.drgbchoices.no_electrodes
                            ii_oth=handles.drgbchoices.no_electrodes+1;
                            this_plv=zeros(1,length(handles_out.PLV(this_ii_PLV).time));
                            if ii_ref<ii_oth
                                this_plv(1,:)=handles_out.PLV(this_ii_PLV).plv(:,ii_ref,ii_oth);
                            else
                                this_plv(1,:)=handles_out.PLV(this_ii_PLV).plv(:,ii_oth,ii_ref);
                            end
                            ii_PLVs=ii_PLVs+1;
                            these_PLVs(ii_PLVs,:)=this_plv;
                            
                        end
                        pls=prctile(these_PLVs(:),95);
                        plot([time_dec(1) time_dec(end)],[pls pls],'-k')
                        
                        ylim([0.2 1])
                        xlabel('Time (sec)')
                        ylabel('PLV')
                        title([handles.drgbchoices.bwlabels{bw_ii} ' ' handles.drgbchoices.percent_labels{per_ii}])
                        
                        if (bw_ii==1)&(per_ii==size(handles.drgbchoices.percent_windows,1))
                            text(0.95,-1.5,'S+','Color','r')
                            text(0.9,-1.5,'S-','Color','b')
                            text(0.9,-1.5,'p<0.05','Color','k')
                        end
                        
                    end
                end
            end
            
            suptitle(['Mouse ' num2str(mouse_ii) 'Phase-locking value'])
            
            
            
            %Plot delta phase timecourse between brain regions
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            
            
            hFig = figure(figNo);
            set(hFig, 'units','normalized','position',[.05 .05 .45 .7])
            
            this_sub=0;
            
            for bw_ii=1:length(handles.drgbchoices.lowF)
                
                for per_ii=size(handles.drgbchoices.percent_windows,1):-1:1
                    this_sub=this_sub+1;
                    subplot(length(handles.drgbchoices.lowF),size(handles.drgbchoices.percent_windows,1),this_sub)
                    hold on
                    for eventNo_ii=length(handles.drgbchoices.evTypeNos):-1:1
                        
                        %Get the PLVs
                        for ii=1:ii_PLV
                            if (handles_out.PLV(ii).eventNo_ii==eventNo_ii)&(handles_out.PLV(ii).per_ii==per_ii)&(handles_out.PLV(ii).bw_ii==bw_ii)
                                this_ii_PLV=ii;
                            end
                        end
                        
                        %Get PLVs
                        these_DPs=zeros(length(handles.drgbchoices.reference_electrodes)*length(handles.drgbchoices.other_electrodes),length(handles_out.PLV(this_ii_PLV).time));
                        ii_PLVs=0;
                        for ii_ref=handles.drgbchoices.reference_electrodes
                            for ii_oth=handles.drgbchoices.other_electrodes
                                this_dp=zeros(1,length(handles_out.PLV(this_ii_PLV).time));
                                if ii_ref<ii_oth
                                    this_dp(1,:)=handles_out.PLV(this_ii_PLV).delta_phase(:,ii_ref,ii_oth);
                                else
                                    this_dp(1,:)=handles_out.PLV(this_ii_PLV).delta_phase(:,ii_oth,ii_ref);
                                end
                                ii_PLVs=ii_PLVs+1;
                                these_DPs(ii_PLVs,:)=this_dp;
                            end
                        end
                        
                        these_DPs_dec=zeros(size(these_DPs,1),size(these_DPs,2)/no_time_points);
                        time_dec=zeros(1,size(these_DPs,2)/no_time_points);
                        time=handles_out.PLV(this_ii_PLV).time;
                        
                        for elec_pair=1:size(these_DPs,1)
                            for ii_dec_tp=1:size(these_DPs,2)/no_time_points
                                these_DPs_dec(elec_pair,ii_dec_tp)=circ_mean(these_DPs(elec_pair,1+(ii_dec_tp-1)*no_time_points:ii_dec_tp*no_time_points)');
                                if elec_pair==1
                                    time_dec(1,ii_dec_tp)=mean(time(1+(ii_dec_tp-1)*no_time_points:ii_dec_tp*no_time_points));
                                end
                            end
                        end
                        
                        mean_dp=circ_mean(these_DPs_dec);
                        
                        
                        if eventNo_ii==2
                            plot(time_dec,mean_dp, 'b');
                        else
                            plot(time_dec,mean_dp,'r');
                        end
                        
                        %                 ylim([0.2 1])
                        xlabel('Time (sec)')
                        ylabel('Delta phase')
                        title([handles.drgbchoices.bwlabels{bw_ii} ' ' handles.drgbchoices.percent_labels{per_ii}])
                        
                        %                 if (bw_ii==1)&(per_ii==size(handles.drgbchoices.percent_windows,1))
                        %                     text(0.95,-1.5,'S+','Color','r')
                        %                     text(0.9,-1.5,'S-','Color','b')
                        %                     text(0.9,-1.5,'p<0.05','Color','k')
                        %                 end
                        
                    end
                end
            end
            
            suptitle(['Mouse ' num2str(mouse_ii) 'Delta phase'])
            
        end
        
    end
    
end

pfft=1;







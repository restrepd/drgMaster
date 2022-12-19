function drgRunBatchMultiDayLFPpar

%drgRunBatchMultiDayLFPpar performs batch analysis of multi day LFP
%
% handles.drgbchoices.analyses chooses the analysis performed
%
% 1 wavelet LFP power
%
% 2 FFT LFP power

close all 
clear all

f_label{1}='Delta';
f_label{2}='Theta';
f_label{3}='Beta';
f_label{4}='Gamma_low';
f_label{5}='Gamma_high';

tic

figNo=0;

first_file=1;

[choiceFileName,choiceBatchPathName] = uigetfile({'drgbChoices*.m'},'Select the .m file with all the choices for analysis');
fprintf(1, ['\ndrgRunBatchMultiDayLFPpar run for ' choiceFileName '\n\n']);

tempDirName=['temp' choiceFileName(12:end-2)];

addpath(choiceBatchPathName)
eval(['handles=' choiceFileName(1:end-2) ';'])
handles.choiceFileName=choiceFileName;
handles.choiceBatchPathName=choiceBatchPathName;

new_no_files=handles.drgbchoices.no_files;
choicePathName=handles.drgbchoices.PathName;
choiceFileName=handles.drgbchoices.FileName;
%Very, very important!
% handles.evTypeNo=handles.drgbchoices.referenceEvent;



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
        fprintf(1, ['Program will be terminated because file No %d does not exist\n'],filNum);
        all_files_present=0;
    end

    if exist( [jtPathName jtFileName(10:end-4) '.edf'])==0
        fprintf(1, ['Program will be terminated because edf file for file No %d does not exist\n'],filNum);
        all_files_present=0;
    end

end


if all_files_present==1

    no_files=handles.drgbchoices.no_files;

    %parfor filNum=first_file:no_files
    for filNum=first_file:no_files
        %         try

        file_no=filNum
        handlespf=struct();
        handlespf=handles;

        this_jt=handlespf.drgbchoices.FileName{filNum};


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

        


        handlespf.drg.drta_p.fullName=[jtPathName jtFileName(10:end-4) '.edf'];

        if handlespf.upload_edf==1
            load([handlespf.drg.drta_p.fullName(1:end-4) '_edf.mat'])
            handlespf.data_out=data_out;
            clear data_out
        end
        
        einf=edfinfo(handlespf.drg.drta_p.fullName);


        %Set the last trial to the last trial in the session
        %handlespf.lastTrialNo=handlespf.drg.session(handlespf.sessionNo).events(2).noTimes;
        handlespf.lastTrialNo=handlespf.drg.session(1).noTrials;

        lastTr=handlespf.lastTrialNo;

        handlespf.burstLowF=min(handles.drgbchoices.lowF);
        handlespf.burstHighF=max(handles.drgbchoices.highF);

        %         handlespf.evTypeNo=1; %Hard coded for continuous
        %         handlespf.subtractRef=0;
        %         handlespf.time_start=-2.8;
        %         handlespf.time_end=5.8;

        edfmat_filename=[handlespf.drg.drta_p.fullName(1:end-4) '_edf.mat'];

        handlespf.ptr_file=matfile(edfmat_filename);


        no_columns=handlespf.ptr_file.no_columns; %hard wired

        mean_dB_power=zeros(length(handlespf.drgbchoices.lowF),lastTr,no_columns);

        do_par=0;
        tic
        process_start_toc=toc;

        %Note: parallel does not work here because having multiple workers call the file
        %pointer slows things down!!!
%         if do_par==1
%             gcp;
%             for chNo=1:no_columns
%                 handlespf.burstLFPNo=chNo;
%                 handlespf.peakLFPNo=chNo;
%                 trialNo=1;
%                 at_end=0;
% 
%                 par_out=[];
%                 for ii=1:ceil(lastTr/100)
%                     par_out(ii).mean_dB_power=zeros(length(handlespf.drgbchoices.lowF),100);
%                 end
% 
%                 %Process 100 trials at a time
%                 %for ii=1:ceil(lastTr/100)
%                 tic
%                 parfor ii=1:ceil(lastTr/100)
%                     tic
%          
%                     handlespff=struct();
%                     handlespff=handlespf;
% 
%                     trialNo=(ii-1)*100+1;
%                     handlespff.lastTrialNo=trialNo+99;
%                     if handlespff.lastTrialNo>lastTr
%                         handlespff.lastTrialNo=lastTr;
%                     end
% 
%                     switch handles.drgbchoices.analyses
%                         case 1
%                             [t_apt,freq,all_Power,all_Power_ref, all_Power_timecourse, this_trialNo, perCorr_pertr, which_event, is_not_moving]=drgGetLFPwavePowerForThisEvTypeNo(handlespff);
%                         case 2
%                             [t,f,all_Power,all_Power_ref, all_Power_timecourse, this_trialNo]=drgGetLFPPowerForThisEvTypeNo(handlespff);
%                     end
% 
%                     for trNo=1:100
%                         frac_not_moving=(sum(is_not_moving(trNo,:))/length(is_not_moving(trNo,:)));
%                         if frac_not_moving>handlespff.drgbchoices.min_frac
%                             log_P_timecourse=zeros(length(freq),length(t_apt));
%                             log_P_timecourse(:,:)=mean(10*log10(all_Power_timecourse),1);
% 
%                             this_mean_dbWB=zeros(1,length(freq));
%                             this_mean_dbWB(1,:)=mean(log_P_timecourse',1);
%                             for ii_freq=1:length(handlespf.drgbchoices.lowF)
%                                 min_f=handlespf.drgbchoices.lowF(ii_freq);
%                                 max_f=handlespf.drgbchoices.highF(ii_freq);
%                                 par_out(ii).mean_dB_power(ii_freq,trNo)=mean(this_mean_dbWB((freq>=min_f)&(freq<=max_f)));
%                             end
%                         else
%                             for ii_freq=1:length(handlespf.drgbchoices.lowF)
%                                 par_out(ii).mean_dB_power(ii_freq,trNo)=NaN;
%                             end
%                         end
%                     end
%                     fprintf(1, ['Processed 100 trial chunk no %d, channel %d\n'], ii,chNo);
%                     fprintf(1, 'Total processing time for this chunk %d minutes\n',(toc)/(60));
%                 end
%                 fprintf(1, 'Total processing time %d hours\n',toc/(60*60));
%                 %Read out the parfor output
%                 for ii=1:ceil(lastTr/100)
%                     this_mean_dB_power=zeros(length(handlespf.drgbchoices.lowF),100);
%                     mean_dB_power(:,(ii-1)*100+1:(ii-1)*100+100,chNo)=par_out(ii).mean_dB_power(:,:);
%                 end
% 
%                 pffft=1;
%             end
%         else
            for chNo=1:no_columns

                trialNo=1;
                at_end=0;
 
                handlespf.burstLFPNo=chNo;
                handlespf.peakLFPNo=chNo;

                %Process 100 trials at a time
                %for ii=1:ceil(lastTr/100)
                for ii=1:ceil(lastTr/100)
                    start_toc=toc;

                    trialNo=(ii-1)*100+1;
                    handlespf.trialNo=trialNo;
                    handlespf.lastTrialNo=trialNo+99;
                    if handlespf.lastTrialNo>lastTr
                        handlespf.lastTrialNo=lastTr;
                    end

                    switch handles.drgbchoices.analyses
                        case 1
                            [t_apt,freq,all_Power,all_Power_ref, all_Power_timecourse, this_trialNo, perCorr_pertr, which_event, is_not_moving]=drgGetLFPwavePowerForThisEvTypeNo(handlespf);
                        case 2
                            [t,freq,all_Power,all_Power_ref, all_Power_timecourse, this_trialNo, perCorr_pertr, which_event, is_not_moving]=drgGetLFPPowerForThisEvTypeNo(handlespf);
                            t_apt=t(1)+[1:size(all_Power_timecourse,3)]*(t(2)-t(1));
                    end

                    for trNo=1:length(this_trialNo)
                        frac_not_moving=(sum(is_not_moving(trNo,:))/length(is_not_moving(trNo,:)));
                        if frac_not_moving>handlespf.drgbchoices.min_frac
                            log_P_timecourse=zeros(length(freq),length(t_apt));
                            log_P_timecourse(:,:)=mean(10*log10(all_Power_timecourse),1);

                            this_mean_dbWB=zeros(1,length(freq));
                            this_mean_dbWB(1,:)=mean(log_P_timecourse',1);
                            for ii_freq=1:length(handlespf.drgbchoices.lowF)
                                min_f=handlespf.drgbchoices.lowF(ii_freq);
                                max_f=handlespf.drgbchoices.highF(ii_freq);
                                mean_dB_power(ii_freq,trNo+trialNo,chNo)=mean(this_mean_dbWB((freq>=min_f)&(freq<=max_f)));
                            end
                        else
                            for ii_freq=1:length(handlespf.drgbchoices.lowF)
                                mean_dB_power(ii_freq,trNo+trialNo,chNo)=NaN;
                            end
                        end
                    end
                    fprintf(1, ['Processed 100 trial chunk no %d, channel %d\n'], ii,chNo);
                    fprintf(1, 'Total processing time for this chunk %d minutes\n',(toc-start_toc)/(60));
                end

                pffft=1;
            end
%         end

        fprintf(1, 'Total overall processing %d hours\n',(toc-process_start_toc)/(60*60));
        
        %Save the data
        save([handlespf.drg.drg_directory handlespf.drg.jt_times_file(10:end-4) handlespf.drgbchoices.save_suffix],'mean_dB_power','handlespf','-v7.3')
        
        
        %Plot the data
        figNo=0;

        
        
        t_hours=handlespf.drgbchoices.start_hour+[1:size(mean_dB_power,2)]*(9/(60*60));

        %Plot dB per bandwidth
        for chNo=1:4
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end

            hFig = figure(figNo);
            set(hFig, 'units','normalized','position',[.07 .05 .45 .75])
            sgtitle(['Power for ' einf.SignalLabels{chNo}])
            for ii_plot=1:5
                subplot(5,1,ii_plot)
                this_mean_dB_power=zeros(1,size(mean_dB_power,2));
                this_mean_dB_power(1,:)=mean_dB_power(ii_plot,:,chNo);
                plot(t_hours,this_mean_dB_power)

                min_dB=prctile(this_mean_dB_power,1)-0.1*(prctile(this_mean_dB_power,99)-prctile(this_mean_dB_power,1));
                max_dB=prctile(this_mean_dB_power,99)+0.1*(prctile(this_mean_dB_power,99)-prctile(this_mean_dB_power,1));
                ylim([min_dB max_dB])
                title(f_label{ii_plot})
                 ylabel('dB')
                xlabel('time (hrs)')


                t=handlespf.drgbchoices.light_on;
              
                %This only works if the start is during the day
                hold on
                while t< t_hours(end)
                    
                    if (rem(t,24)<handlespf.drgbchoices.light_on)||(rem(t,24)>=handlespf.drgbchoices.light_off)
                        %This is the night
                        
                        delta_t_night_end=24-handlespf.drgbchoices.light_off+handlespf.drgbchoices.light_on;
                        if t+delta_t_night_end<t_hours(end)
                            plot([t t+delta_t_night_end], [min_dB+0.1*(prctile(this_mean_dB_power,99)-prctile(this_mean_dB_power,1))...
                                min_dB+0.1*(prctile(this_mean_dB_power,99)-prctile(this_mean_dB_power,1)) ], '-k', 'LineWidth', 3)
                        else
                            plot([t t_hours(end)], [min_dB+0.1*(prctile(this_mean_dB_power,99)-prctile(this_mean_dB_power,1))...
                                min_dB+0.1*(prctile(this_mean_dB_power,99)-prctile(this_mean_dB_power,1)) ], '-k', 'LineWidth', 3)
                        end
                        t=t+delta_t_night_end;

                    else
                        %This is the day
                        t_last=t;
                        t=t+(handlespf.drgbchoices.light_off-rem(t,24));


                    end
                    
                end

            end
        end

        %Subtract delta
        %Plot dB per bandwidth
        for chNo=1:4
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end

            hFig = figure(figNo);
            set(hFig, 'units','normalized','position',[.07 .05 .45 .75])
            sgtitle(['Delta power subtracted for ' einf.SignalLabels{chNo}])
            for ii_plot=2:5
                subplot(4,1,ii_plot-1)
                this_mean_dB_power=zeros(1,size(mean_dB_power,2));
                this_mean_dB_power(1,:)=mean_dB_power(ii_plot,:,chNo)-mean_dB_power(1,:,chNo);
                plot(t_hours,this_mean_dB_power)

                min_dB=prctile(this_mean_dB_power,1)-0.1*(prctile(this_mean_dB_power,99)-prctile(this_mean_dB_power,1));
                max_dB=prctile(this_mean_dB_power,99)+0.1*(prctile(this_mean_dB_power,99)-prctile(this_mean_dB_power,1));
                ylim([min_dB max_dB])
                title(f_label{ii_plot})
                ylabel('delta dB')
                xlabel('time (hrs)')

                t=handlespf.drgbchoices.light_on;
              
                %This only works if the start is during the day
                hold on
                while t< t_hours(end)
                    
                    if (rem(t,24)<handlespf.drgbchoices.light_on)||(rem(t,24)>=handlespf.drgbchoices.light_off)
                        %This is the night
                        
                        delta_t_night_end=24-handlespf.drgbchoices.light_off+handlespf.drgbchoices.light_on;
                        if t+delta_t_night_end<t_hours(end)
                            plot([t t+delta_t_night_end], [min_dB+0.1*(prctile(this_mean_dB_power,99)-prctile(this_mean_dB_power,1))...
                                min_dB+0.1*(prctile(this_mean_dB_power,99)-prctile(this_mean_dB_power,1)) ], '-k', 'LineWidth', 3)
                        else
                            plot([t t_hours(end)], [min_dB+0.1*(prctile(this_mean_dB_power,99)-prctile(this_mean_dB_power,1))...
                                min_dB+0.1*(prctile(this_mean_dB_power,99)-prctile(this_mean_dB_power,1)) ], '-k', 'LineWidth', 3)
                        end
                        t=t+delta_t_night_end;

                    else
                        %This is the day
                        t_last=t;
                        t=t+(handlespf.drgbchoices.light_off-rem(t,24));


                    end
                    
                end

            end
        end
    end

end







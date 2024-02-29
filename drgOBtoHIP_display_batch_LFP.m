function drgOBtoHIP_display_batch_LFP(choiceBatchPathName,choiceFileName)
%Note: fitcnet will not work in Matlab versions earlier than 2021a


first_file=1;

if nargin==0
    [choiceFileName,choiceBatchPathName] = uigetfile({'drgbChoicesOBtoHIP*.m'},'Select the .m file with all the choices for analysis');
end

fprintf(1, ['\ndrgOBtoHIP_batch_LFP run for ' choiceFileName '\n\n']);

addpath(choiceBatchPathName)
eval(['handles=' choiceFileName(1:end-2) ';'])
handles.choiceFileName=choiceFileName;
handles.choiceBatchPathName=choiceBatchPathName;

new_no_files=handles.no_files;


%Parallel batch processing for each file
all_files_present=1;
for filNum=first_file:handles.no_files
     
    
    %Make sure that all the files exist
    jtFileName=handles.FileName{filNum};
    if iscell(handles.PathName)
        jtPathName=handles.PathName{filNum};
    else
        jtPathName=handles.PathName;
    end
     
    if exist([jtPathName jtFileName])==0
        fprintf(1, ['Program will be terminated because file No %d, ' jtPathName jtFileName ' does not exist\n'],filNum);
        all_files_present=0;
    end
    
end



tic

if all_files_present==1
    
    
    %Show spectrograms per mouse per electrode
    these_mice=unique(handles.mouse_no);
    figNo=0;
    for mouseNo=1:length(these_mice)
        %Find the files for this mouse
        these_files=[];
        for fileNo=1:handles.no_files
            if handles.mouse_no(fileNo)==mouseNo
                %Does this file exist?
                PathName_out=handles.PathName{fileNo};
                this_filename=handles.FileName{fileNo};
                FileName_out=['OBtoHIPLFP_' this_filename(11:end)];
                if exist([PathName_out FileName_out])~=0
                    these_files=[these_files fileNo];
                end
            end
        end

        for ii_electrode=1:length(handles.peakLFPNo)
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            hFig2 = figure(figNo);
            set(hFig2, 'units','normalized','position',[.1 .1 .7 .7])
            hold on

            % suptitle(['Spectrograms for ' handles.electrode_label{handles.peakLFPNo(ii_electrode)} handles.MouseName{mouseNo}])
            %Get max and min
            all_log_P_timecourse=[];
            for fileNo=these_files
                %Plot the timecourse
                PathName_out=handles.PathName{fileNo};
                this_filename=handles.FileName{fileNo};
                FileName_out=['OBtoHIPLFP_' this_filename(11:end)];
                load([PathName_out FileName_out])

                
                log_P_timecourse=handles_out.electorde(1).handles_out.log_P_timecourse;
                all_log_P_timecourse=[all_log_P_timecourse log_P_timecourse(:)'];
                
            end
            maxLogPper=prctile(all_log_P_timecourse(:),99);
        minLogPper=prctile(all_log_P_timecourse(:),1);
        sessionNo=0;
            for fileNo=these_files
                %Plot the timecourse
                sessionNo=sessionNo+1;
                PathName_out=handles.PathName{fileNo};
                this_filename=handles.FileName{fileNo};
                FileName_out=['OBtoHIPLFP_' this_filename(11:end)];
                load([PathName_out FileName_out])

                time=handles_out.electorde(1).handles_out.time;
                log_P_timecourse=handles_out.electorde(1).handles_out.log_P_timecourse;
                f=handles_out.electorde(1).handles_out.f;

                subplot(length(these_files),1,sessionNo)

                drg_pcolor(repmat(time,length(f),1)',repmat(f,length(time),1),log_P_timecourse')

                colormap fire
                shading interp
                caxis([minLogPper maxLogPper]);

                hold on
                %     plot([time(handles.ii_laser_start) time(handles.ii_laser_start)],[f(1) f(end)],'-k','LineWidth',2)
                %     plot([time(handles.ii_laser_end) time(handles.ii_laser_end)],[f(1) f(end)],'-k','LineWidth',2)

                if fileNo==length(these_files)
                xlabel('Time (min)')
                end
                ylabel('Frequency (Hz)');
                title([handles.MouseName{mouseNo} ' ' handles.electrode_label{handles.peakLFPNo(ii_electrode)} ' session No ' num2str(sessionNo)])
            end
        end
    end

     
end

pffft=1;





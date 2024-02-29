function drgOBtoHIP_batch_LFP(choiceBatchPathName,choiceFileName)
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
    
    
    %Process each file separately
    at_end=0;
    fileNo=0;
    while at_end==0

        if fileNo+1>handles.no_files
            at_end=1;
        else
            %Process the next session
            first_toc=toc;

            handles_in=[];

            fileNo=fileNo+1;

            ii_f=1;
            this_session=handles.session_no(fileNo);
            handles_in.jtPathNames{ii_f}=handles.PathName{fileNo};
            handles_in.jtFileNames{ii_f}=handles.FileName{fileNo};

            handles.PathName_out=handles.PathName{fileNo};
            this_filename=handles.FileName{fileNo};
            handles.FileName_out=['OBtoHIPLFP_' this_filename(11:end)];

            handles_in.ii_laser_start=handles.ii_laser_start(fileNo);
            handles_in.ii_laser_end=handles.ii_laser_end(fileNo);

            found_session=0;
            this_session=handles.session_no(fileNo);
            while found_session==0
                if fileNo==handles.no_files
                    found_session=1;
                else
                    if handles.session_no(fileNo+1)==this_session
                        ii_f=ii_f+1;
                        fileNo=fileNo+1;
                        handles_in.jtPathNames{ii_f}=handles.PathName{fileNo};
                        handles_in.jtFileNames{ii_f}=handles.FileName_jt{fileNo};
                    else
                        found_session=1;
                    end
                end
            end
            
            handles_out=[];
            handles_in.electrode_label=handles.electrode_label;

            handles_in.window=handles.window;
            handles_in.noverlap=handles.noverlap;

            handles_in.burstLowF=handles.burstLowF;
            handles_in.burstHighF=handles.burstHighF;
            handles_in.showData=handles.showData;


            for ii_electrode=1:length(handles.peakLFPNo)
                handles_in.peakLFPNo=handles.peakLFPNo(ii_electrode);
                handles_out.electorde(ii_electrode).handles_out=drgOBtoHIP_LFPpower(handles_in);
            end


            fprintf(1,'\n\nProcessing time for file No %d is %d hours\n',fileNo,(toc-first_toc)/(60*60));

            %Save output file
            handles_out.handles=handles;
            save([handles.PathName_out handles.FileName_out],'handles_out','-v7.3')
        end
    end
    
   
    
    fprintf(1, 'Total processing time %d hours\n',toc/(60*60));
     
end







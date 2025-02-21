function drgOBtoHIP_batch_tPRP(choiceBatchPathName,choiceFileName)
%Note: fitcnet will not work in Matlab versions earlier than 2021a




if nargin==0
    [choiceFileName,choiceBatchPathName] = uigetfile({'drgbChoicesOBtoHIPtPRP*.m'},'Select the .m file with all the choices for analysis');
end

fprintf(1, ['\ndrgOBtoHIP_batch_LFP run for ' choiceFileName '\n\n']);

addpath(choiceBatchPathName)
eval(['handles=' choiceFileName(1:end-2) ';'])
handles.choiceFileName=choiceFileName;
handles.choiceBatchPathName=choiceBatchPathName;

if handles.is_sphgpu==1
    addpath('/home/restrepd/Documents/MATLAB/drgMaster')
    addpath('/home/restrepd/Documents/MATLAB/m new/CircStat2012a')
end

new_no_files=handles.no_files;


%Parallel batch processing for each file
all_files_present=1;
first_file=handles.first_file;
for filNum=first_file:handles.no_files
     
    
    %Make sure that all the files exist
    jtFileName=handles.FileName{filNum};
    if iscell(handles.PathName)
        jtPathName=handles.PathName{filNum};
    else
        jtPathName=handles.PathName;
    end
     
    if exist([jtPathName jtFileName])==0
        fprintf(1, ['Program will be terminated because file No %d, ' jtFileName ' does not exist\n\n'],filNum);
        all_files_present=0;
    end
    
end



tic

if all_files_present==1


    %Process each file separately
    for fileNo=first_file:handles.no_files
        
        handles.PathName_out=handles.PathName{fileNo};
        this_filename=handles.FileName{fileNo};

        handles.FileName_out=['OBtoHIPtPRP_' this_filename(11:end)];

        if (exist([handles.PathName_out handles.FileName_out])==2)&(handles.overwrite_files==0)
            fprintf(1,'\n\nOutput file for file No %d already exists and will not be overwritten\n\n',fileNo);
        else
            %Process the next file
            first_toc=toc;

            handles_in=[];
            handles_out=[];

            handles_in.jtPathNames{1}=handles.PathName{fileNo};
            handles_in.jtFileNames{1}=handles.FileName{fileNo};

            handles_in.electrode_label=handles.electrode_label;

            handles_in.window=handles.window;
            handles_in.noverlap=handles.noverlap;

            handles_in.burstLowF=handles.burstLowF;
            handles_in.burstHighF=handles.burstHighF;
            handles_in.showData=handles.showData;
            handles_in.peakLowF=handles.peakLowF;
            handles_in.peakHighF=handles.peakHighF;

            handles_in.n_phase_bins=handles.n_phase_bins;

            handles_in.which_method=handles.which_method;
            handles_in.save_drgb=handles.save_drgb;
            handles_in.use_peakAngle=handles.use_peakAngle;

            handles_in.window=handles.window; %This is the FFT window in sec  %Tort is 1 sec, old DR 0.37
            handles_in.noverlap=handles.noverlap;
            
            for ii_electrode=1:length(handles.peakLFPNo)
                fprintf(1,'\n\nProcessing file No %d, electrode %d\n\n',fileNo,handles.peakLFPNo(ii_electrode));
                handles_in.peakLFPNo=handles.peakLFPNo(ii_electrode);
                handles_out.electrode(ii_electrode).handles_out=drgOBtoHIP_tPRP(handles_in);
                fprintf(1,'\n\nProcessed file No %d, electrode %d\n\n',fileNo,handles.peakLFPNo(ii_electrode));
            end

            fprintf(1,'\n\nProcessing time for file No %d is %d hours\n\n',fileNo,(toc-first_toc)/(60*60));

            %Save output file
            handles_out.handles=handles;
            save([handles.PathName_out handles.FileName_out],'handles_out','-v7.3')
        end

    end

    fprintf(1, 'Total processing time %d hours\n',toc/(60*60));
     
end







function drgChR2_batch_LFPpower(choiceBatchPathName,choiceFileName)
%Note: fitcnet will not work in Matlab versions earlier than 2021a

if nargin==0
    [choiceFileName,choiceBatchPathName] = uigetfile({'drgChR2_choices*.m'},'Select the .m file with all the choices for analysis');
end

fprintf(1, ['\nddrgChR2_batch_LFPpower run for ' choiceFileName '\n\n']);

tempDirName=['temp' choiceFileName(12:end-2)];

addpath(choiceBatchPathName)
eval(['handles=' choiceFileName(1:end-2) ';'])
handles.choiceFileName=choiceFileName;
handles.choiceBatchPathName=choiceBatchPathName;

%Are all files there?
all_files_present=1;
first_file=handles.first_file;
for filNum=first_file:handles.no_files
     
    
    %Make sure that all the files exist
    FileName_in=handles.FileName_in{filNum};
    if iscell(handles.PathName_in)
        PathName_in=handles.PathName_in{filNum};
    else
        PathName_in=handles.PathName_in;
    end
      
    if exist([PathName_in FileName_in])==0
        fprintf(1, ['Program will be terminated because file No %d, ' PathName_in FileName_in ' does not exist\n'],filNum);
        all_files_present=0;
    end
    
end




%process each file for each electrode
if all_files_present==1

    
    
    this_handles_choices.burstLowF=handles.burstLowF;
    this_handles_choices.burstHighF=handles.burstHighF;
    
    


    %Process each file separately
    tic

    for fileNo=first_file:length(handles.FileName_in)
%         this_handles_choices.peakLFPNo=handles.file(fileNo).peakLFPNo
        this_handles_choices.jtFileName=handles.FileName_in{fileNo};
        this_handles_choices.jtPathName=handles.PathName_in{fileNo};
        handles_per_file=[];

        this_ii=0;
        for peakLFPNo_ii=1:length(handles.file(fileNo).peakLFPNo)
            this_ii=this_ii+1;
            this_handles_choices.electrode_label=handles.electrode_labels{this_ii};
            this_handles_choices.peakLFPNo=handles.file(fileNo).peakLFPNo(peakLFPNo_ii);
            start_toc=toc;
            handles_per_file.peakLFPNo(peakLFPNo_ii).handles_out=drgOBChr2LFPpower(this_handles_choices);
            fprintf(1, ['Data processed for file number %d took %d minutes\n'],fileNo,(toc-start_toc)/60);
        end
    

     %Save output file
 
        save([handles.PathName_out FileName_in(1:end-4) handles.suffix_out],'handles_per_file','handles','-v7.3')
    end
    fprintf(1, 'Total processing time %d hours\n',toc/(60*60));

end







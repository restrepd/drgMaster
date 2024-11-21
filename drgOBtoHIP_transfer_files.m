function drgOBtoHIP_transfer_files(choiceBatchPathName,choiceFileName)
%Note: fitcnet will not work in Matlab versions earlier than 2021a


if nargin==0
    [choiceFileName,choiceBatchPathName] = uigetfile({'drgbChoicesOBtoHIP*.m'},'Select the .m file with all the choices for analysis');
end

fprintf(1, ['\ndrgOBtoHIP_batch_LFP run for ' choiceFileName '\n\n']);

addpath(choiceBatchPathName)
eval(['handles=' choiceFileName(1:end-2) ';'])
handles.choiceFileName=choiceFileName;
handles.choiceBatchPathName=choiceBatchPathName;

out_PathName='/data2/5xFADvsWT_1_hour_treatmentVsnone/OBtoHIP_out/';

try
    mkdir(out_PathName)
catch
end

%Transfer output files
for filNum=1:handles.no_files
     
    if iscell(handles.PathName)
        jtPathName=handles.PathName{filNum};
    else
        jtPathName=handles.PathName;
    end

    this_filename=handles.FileName{filNum};

    FileName_out=['OBtoHIPLFP_' this_filename(11:end)];

    copyfile([jtPathName FileName_out],[out_PathName FileName_out])

    fprintf(1, ['Transfered LFP file No %d, ' FileName_out '\n'],filNum);

    FileName_out=['OBtoHIPtPRP_' this_filename(11:end)];

    copyfile([jtPathName FileName_out],[out_PathName FileName_out])

    fprintf(1, ['Transfered LFP file No %d, ' FileName_out '\n'],filNum);

end



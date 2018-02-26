function drgRunBatchSpikes

%Ask user for the m file that contains information on what the user wants the analysis to be
%This file has all the information on what the user wants done, which files
%to process, what groups they fall into, etc
%
% An example of this file: drgbChoicesDanielPrelim
%
%


[choiceFileName,choiceBatchPathName] = uigetfile({'drgbChoices*.m'},'Select the .m file with all the choices for analysis');
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
first_file=handles.drgb.first_file;

if first_file==1
    handles.drgb.unit_no=0;
else
    load([handles.drgb.outPathName handles.drgb.outFileName])
    handles.drgb=handles_drgb.drgb;
    %The user may add new files
    handles.drgbchoices.no_files=new_no_files;
    handles.drgbchoices.PathName=choicePathName;
    handles.drgbchoices.FileName=choiceFileName;
end

test_batch=handles.drgbchoices.test_batch;

%Do batch processing for each file
for filNum=first_file:handles.drgbchoices.no_files
    
    file_no=filNum
    
    if test_batch==1
        if handles.drgbchoices.group_no(filNum)==1
            handles.data_vs_simulate=1;
        else
            handles.data_vs_simulate=4;
        end
    end
    
    %read the jt_times files
    handles.jtfullName=[handles.drgbchoices.PathName{filNum},handles.drgbchoices.FileName{filNum}];
    handles.jtFileName=handles.drgbchoices.FileName{filNum};
    handles.jtPathName=handles.drgbchoices.PathName{filNum};
    
    drgRead_jt_times(handles.jtPathName,handles.jtFileName);
    FileName=[handles.jtFileName(10:end-4) '_drg.mat'];
    handles.fullName=[handles.jtPathName,FileName];
    handles.FileName=FileName;
    handles.PathName=handles.jtPathName;
    
    load(handles.fullName);
    handles.drg=drg;
    
    if handles.read_entire_file==1
        handles=drgReadAllDraOrDg(handles);
    end
    
    switch handles.drg.session(handles.sessionNo).draq_p.dgordra
        case 1
        case 2
            handles.drg.drta_p.fullName=[handles.jtPathName handles.jtFileName(10:end-4) '.dg'];
        case 3
            handles.drg.drta_p.fullName=[handles.jtPathName handles.jtFileName(10:end-4) '.rhd'];
    end
    
    
    
    %Set the last trial to the last trial in the session
    handles.lastTrialNo=handles.drg.session(handles.sessionNo).events(2).noTimes;
    
    %Save information for this file
    handles.drgb.filNum=filNum;
    handles.drgb.file(filNum).FileName=handles.FileName;
    handles.drgb.file(filNum).PathName=handles.PathName;
    
    
    handles.drgb.file(filNum).drg=handles.drg;
    
    %Run the analysis for a time range encompassing all windows
    
        
        handles.time_start=handles.drgbchoices.minStart-handles.time_pad;
        handles.time_end=handles.drgbchoices.maxEnd+handles.time_pad; 
        
        
        %Now run the analysis for each unit
        for unitNo=1:handles.drg.session(handles.sessionNo).noUnits
            handles.unitNo=unitNo;
            handles.drgb.unit_no=handles.drgb.unit_no+1;
            handles.drgb.unit(handles.drgb.unit_no).SingleUnit=handles.drg.unit(unitNo).SingleUnit;
            handles.drgb.unit(handles.drgb.unit_no).SingleUnit_Fee=handles.drg.unit(unitNo).SingleUnit_Fee;
            handles.drgb.unit(handles.drgb.unit_no).basalFR=handles.drg.unit(unitNo).basalFR;
            handles.drgb.unit(handles.drgb.unit_no).perViol=handles.drg.unit(unitNo).perViol;
            handles.drgb.unit(handles.drgb.unit_no).groupNo=handles.drgbchoices.group_no(filNum);
            handles=drgPlotPSTHFct(handles);
        end
        
    
    
    %Save output file
    handles_drgb=handles;
    if isfield(handles,'data_dg')
        handles_drgb=rmfield(handles_drgb,'data_dg');
    end
    save([handles.drgb.outPathName handles.drgb.outFileName],'handles_drgb','-v7.3')
    
end





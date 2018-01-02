function drgRunBatchBehavior

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
    handles.drgb.lfpevpair_no=0;
    handles.drgb.lfp_per_exp_no=0;
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
            handles.data_vs_simulate=5;
        else
            handles.data_vs_simulate=6;
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
    
    [handles.drgb.file(filNum).perCorr, handles.drgb.file(filNum).encoding_trials, handles.drgb.file(filNum).retrieval_trials, encoding_this_evTypeNo,retrieval_this_evTypeNo]=drgFindEncRetr(handles);


    
%     handles.drgb.file(filNum).drg=handles.drg;
    
    
%     %Save output file
%     handles_drgb=handles;
%     if isfield(handles,'data_dg')
%         handles_drgb=rmfield(handles_drgb,'data_dg');
%     end
%     save([handles.drgb.outPathName handles.drgb.outFileName],'handles_drgb','-v7.3')
    
end
 
%Plot percent correct
try
    close 1
catch
end

hFig1 = figure(1);
set(hFig1, 'units','normalized','position',[.02 .02 .95 .95])

max_session=max(handles.drgbchoices.session_no);
max_mouse=max(handles.drgbchoices.mouse_no);

for filNum=first_file:handles.drgbchoices.no_files
    subplot(max_mouse,max_session,max_session*(handles.drgbchoices.mouse_no(filNum)-1)+handles.drgbchoices.session_no(filNum))
    % subplot(3,1,1)
    trials=1:length(handles.drgb.file(filNum).perCorr);
    
    %Plot in different colors
    plot(trials,handles.drgb.file(filNum).perCorr,'o','MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7])
    hold on
    plot(trials(handles.drgb.file(filNum).encoding_trials),handles.drgb.file(filNum).perCorr(handles.drgb.file(filNum).encoding_trials),'ob')
    plot(trials(handles.drgb.file(filNum).retrieval_trials),handles.drgb.file(filNum).perCorr(handles.drgb.file(filNum).retrieval_trials),'or')
    
    ylim([0 110]);
    title(handles.drgbchoices.group_no_names{handles.drgbchoices.group_no(filNum)})
end

title_str=inputdlg('Enter title');

annotation('textbox', [0 0.9 1 0.1], ...
    'String', title_str{1}, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
pfft=1;






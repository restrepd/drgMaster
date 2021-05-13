function drgRunBatchLFPCrossCorrpar

%drgRunBatchLFPCorssCorrpar performs batch analysis of trial per trial lag 
%following the procedure of Adhikari et al. 2010
%https://doi.org/10.1016/j.jneumeth.2010.06.019
%Figures 4F,G

tic

first_file=1;

[choiceFileName,choiceBatchPathName] = uigetfile({'drgbChoicesCrossCorr*.m'},'Select the .m file with all the choices for analysis');
fprintf(1, ['\ndrgRunBatchLFPCrossCorrpar run for ' choiceFileName '\n\n']);

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
handles_corr=[];
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
    this_jt=handles.drgbchoices.FileName{filNum};
    handles.temp_exist(filNum)=exist([handles.drgbchoices.PathName{filNum} 'Corr__' this_jt(10:end)]);
    
    if handles.temp_exist(filNum)==2
        %If it was processed load the temp result
        load([handles.drgbchoices.PathName{filNum} 'Corr__' this_jt(10:end)])
        handles_corr(filNum).corr_out=handles_corr_out;
    else
        handles_corr(filNum).corr_out=[];
    end
end


if all_files_present==1
    
    
    no_files=handles.drgbchoices.no_files;
    
    for filNum=first_file:no_files
        %             for filNum=first_file:no_files
        %         try
        
        file_no=filNum
        
        
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
            
            switch handles.drg.session(handles.sessionNo).draq_p.dgordra
                case 1
                case 2
                    handles.drg.drta_p.fullName=[jtPathName jtFileName(10:end-4) '.dg'];
                case 3
                    handles.drg.drta_p.fullName=[jtPathName jtFileName(10:end-4) '.rhd'];
            end
            
            handles.time_start=handles.drgbchoices.time_start;
            handles.time_end=handles.drgbchoices.time_end;
            handles.peakLowF=handles.drgbchoices.lowF;
            handles.peakHighF=handles.drgbchoices.highF;
            handles.burstLowF=handles.drgbchoices.lowF;
            handles.burstHighF=handles.drgbchoices.highF;
            handles.jtPathName=jtPathName;
            handles.jtFileName=jtFileName;
            handles.trialNo=1;
            handles.lastTrialNo=handles.drg.session(1).noTrials;
    
            [handles, handles_corr_out]=drgLFPCorrTimecourseAllTets(handles);
            
            handles_corr(filNum).corr_out=handles_corr_out;
            
            if (exist([handles.jtPathName 'Corr_' handles.jtFileName(9:end)])==2)
                fprintf(1, 'Loaded pre-processed cross-correlation for file number: %d\n',filNum);
            else
                fprintf(1, 'Processed cross-correlation for file number: %d\n',filNum);
            end

    end
    
    
end


%Save output file
handles_drgb=handles;
if isfield(handles,'data_dg')
    handles_drgb=rmfield(handles_drgb,'data_dg');
end
save([handles.drgb.outPathName handles.drgb.outFileName],'handles_drgb','-v7.3')

fprintf(1, 'Total processing time %d hours\n',toc/(60*60));

end







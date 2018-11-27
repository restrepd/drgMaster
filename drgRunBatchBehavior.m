function drgRunBatchBehavior

%Ask user for the m file that contains information on what the user wants the analysis to be
%This file has all the information on what the user wants done, which files
%to process, what groups they fall into, etc
%
% An example of this file: drgbChoicesDanielPrelim
%
%
pffft = 1;

close all
clear all

%
% % which_display=1 used to show first and last behavior
% % For Daniel's new Fig. 1
% % drgbChoicesDanielAPEBfirstandlastBeh11618.m
% % drgbChoicesDanielEAPAfirstandlastBeh11418.m
% % drgbChoicesDanielIAMOirstandlastBeh11618.m
% which_display=1;
% trial_window=30;

% which_display=2 is used to show behavior for per and post laser for Fig.
% 6 of Daniel's paper
% drgbChoicesDanielAPEBexperimentalBeh02012018
which_display=1;
trial_window=20;

which_file=1; %1=.m   2=.mat

if which_file==1
    [choiceFileName,choiceBatchPathName] = uigetfile({'drgbChoices*.m'},'Select the .m file with all the choices for analysis');
    addpath(choiceBatchPathName)
    eval(['handles=' choiceFileName(1:end-2) ';'])
    handles.choiceFileName=choiceFileName;
    handles.choiceBatchPathName=choiceBatchPathName;
    
    new_no_files=length(handles.drgbchoices.PathName);
    choicePathName=handles.drgbchoices.PathName;
    choiceFileName=handles.drgbchoices.FileName;
    
    %Do batch processing for each file
    for filNum=1:length(handles.drgbchoices.FileName)
        
        file_no=filNum
        
        this_file=handles.drgbchoices.FileName{filNum};
        
        if strcmp(this_file(1:3),'jt_')
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
            
        else
            %Read dropc .mat file
            handles.dropc_hand=drg_dropc_load([handles.drgbchoices.PathName{filNum},handles.drgbchoices.FileName{filNum}]);
            
            
            %Compute percent correct in a 20 trial window
            sliding_window=20; %Trials for determination of behavioral performance
            
            
            no_trials=length(handles.dropc_hand.dropcData.trialTime);
            score=~(handles.dropc_hand.dropcData.trialScore==(handles.dropc_hand.dropcData.odorType-1));
            
            for ii=1:length(handles.dropc_hand.dropcData.trialTime)-sliding_window+1
                first_time=handles.dropc_hand.dropcData.trialTime(ii);
                last_time=handles.dropc_hand.dropcData.trialTime(ii+sliding_window-1);
                handles.drgb.file(filNum).perCorr(ii+(sliding_window/2))=100*sum(score(ii:ii+sliding_window-1))/sliding_window;
                
                if ii==1
                    handles.drgb.file(filNum).perCorr(ii:(sliding_window/2))=handles.drgb.file(filNum).perCorr(ii+(sliding_window/2));
                end
                if ii==length(handles.dropc_hand.dropcData.trialTime)-sliding_window+1
                    handles.drgb.file(filNum).perCorr(ii+(sliding_window/2)+1:length(handles.dropc_hand.dropcData.trialTime))=handles.drgb.file(filNum).perCorr(ii+(sliding_window/2));
                end
            end
            
            handles.drgb.file(filNum).encoding_trials=handles.drgb.file(filNum).perCorr<=65;
            handles.drgb.file(filNum).retrieval_trials=handles.drgb.file(filNum).perCorr>=80;
            
        end
        
        
    end
    
    %Save the data
    save([handles.choiceBatchPathName handles.choiceFileName(1:end-2) '.mat'],'handles')
    
else
    [matFileName,matBatchPathName] = uigetfile({'drgbChoices*.mat'},'Select the .mat file with all the choices for analysis');
    load([matBatchPathName matFileName])
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

for filNum=1:length(handles.drgbchoices.FileName)
    subplot(max_mouse,max_session,max_session*(handles.drgbchoices.mouse_no(filNum)-1)+handles.drgbchoices.session_no(filNum))
    set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
    % subplot(3,1,1)
    trials=1:length(handles.drgb.file(filNum).perCorr);
    
    %Plot in different colors
    plot(trials,handles.drgb.file(filNum).perCorr,'o','MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7])
    hold on
    plot(trials(handles.drgb.file(filNum).encoding_trials),handles.drgb.file(filNum).perCorr(handles.drgb.file(filNum).encoding_trials),'ob')
    plot(trials(handles.drgb.file(filNum).retrieval_trials),handles.drgb.file(filNum).perCorr(handles.drgb.file(filNum).retrieval_trials),'or')
    
    ylim([0 110]);
    title([handles.drgbchoices.group_no_names{handles.drgbchoices.group_no(filNum)} ':' handles.drgbchoices.epoch_names{handles.drgbchoices.epoch(filNum)}])
    
    remainder = rem((max_session*(handles.drgbchoices.mouse_no(filNum)-1)+handles.drgbchoices.session_no(filNum))-1,max_session);
    if  remainder==0
        ylabel(handles.drgbchoices.MouseName(handles.drgbchoices.mouse_no(filNum)))
    end
end

title_str=inputdlg('Enter title');

annotation('textbox', [0 0.9 1 0.1], ...
    'String', title_str{1}, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')



switch which_display
    case 1
        %Plot the percent correct for the first and last set of trials
        
        first_pc=[];
        first_pc_end=[];
        last_pc=[];
        
        for filNum=1:length(handles.drgbchoices.FileName)
            if length(handles.drgb.file(filNum).perCorr)>=trial_window
                if handles.drgbchoices.session_no(filNum)==1
                    first_pc(handles.drgbchoices.mouse_no(filNum))=mean(handles.drgb.file(filNum).perCorr(1:trial_window));
                    first_pc_end(handles.drgbchoices.mouse_no(filNum))=mean(handles.drgb.file(filNum).perCorr(end-trial_window:end));
                else
                    if handles.drgbchoices.session_no(filNum)==max(handles.drgbchoices.session_no)
                        last_pc(handles.drgbchoices.mouse_no(filNum))=mean(handles.drgb.file(filNum).perCorr(end-trial_window:end));
                    end
                end
            end
        end
        
        %Plot pc of the first 30 trials in first and last 30 trials in last
        figure(2)
        hold on
        for fps=1:max(handles.drgbchoices.mouse_no)
            plot([0 1],[first_pc(fps) last_pc(fps)],'-o', 'Color',[0.7 0.7 0.7])
        end
        
        plot([0 1],[mean(first_pc) mean(last_pc)],'-k','LineWidth', 3)
        CI = bootci(1000, @mean, first_pc);
        plot([0 0],CI,'-b','LineWidth',3)
        plot(0,mean(first_pc),'ob','MarkerSize', 10,'MarkerFace','b')
        CI = bootci(1000, @mean, last_pc);
        plot([1 1],CI,'-r','LineWidth',3)
        plot(1,mean(last_pc),'or','MarkerSize', 10,'MarkerFace','r')
        ylabel('Percent correct')
        ylim([30 110])
        title('Percent correct, first 30 trials, last 30')
        set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
        
        p_perCorr=ranksum(last_pc,first_pc);
        fprintf(1, '\np value for ranksum test for percent correct= %d\n\n',p_perCorr);
        
        
        %Plot one of the sessions
        filNum=6;
        figure(4)
        set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
        % subplot(3,1,1)
        trials=1:length(handles.drgb.file(filNum).perCorr);
        
        %Plot in different colors
        plot(trials,handles.drgb.file(filNum).perCorr,'o','MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7])
        hold on
        plot(trials(handles.drgb.file(filNum).encoding_trials),handles.drgb.file(filNum).perCorr(handles.drgb.file(filNum).encoding_trials),'ob')
        plot(trials(handles.drgb.file(filNum).retrieval_trials),handles.drgb.file(filNum).perCorr(handles.drgb.file(filNum).retrieval_trials),'or')
        
        ylim([30 110]);
        title(['Mouse: ' handles.drgbchoices.MouseName{handles.drgbchoices.mouse_no(filNum)} ', session No: ' num2str(handles.drgbchoices.session_no(filNum))])
        
        set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
        pffft=1
        
    case 2
        %Plot the percent correct for the last 30 trials in pre and post laser
        
        for grNo=1:max(handles.drgbchoices.group_no)
            
            
            group(grNo).first_pc_end=[];
            
            group(grNo).last_pc_start=[];
            
            group(grNo).mouse_numbers=[];
            for filNum=1:length(handles.drgbchoices.FileName)
                if length(handles.drgb.file(filNum).perCorr)>=trial_window
                    if handles.drgbchoices.group_no(filNum)==grNo
                        if handles.drgbchoices.session_no(filNum)==1
                            group(grNo).first_pc_end(handles.drgbchoices.mouse_no(filNum))=mean(handles.drgb.file(filNum).perCorr(end-trial_window:end));
                            group(grNo).mouse_numbers=[group(grNo).mouse_numbers handles.drgbchoices.mouse_no(filNum)];
                        else
                            if handles.drgbchoices.session_no(filNum)==max(handles.drgbchoices.session_no)
                                group(grNo).last_pc_start(handles.drgbchoices.mouse_no(filNum))=mean(handles.drgb.file(filNum).perCorr(1:trial_window));
                            end
                        end
                    end
                end
            end
            
            %Plot pc of the first 30 trials in first and last 30 trials in last
            
            try
                close(grNo+1)
            catch
            end
            figure(grNo+1)
            hold on
            for filNum=1:length(handles.drgbchoices.FileName)
                if length(handles.drgb.file(filNum).perCorr)>=trial_window
                    if handles.drgbchoices.group_no(filNum)==grNo
                        
                        plot([0 1],[group(grNo).first_pc_end(handles.drgbchoices.mouse_no(filNum)) group(grNo).last_pc_start(handles.drgbchoices.mouse_no(filNum))],'-o', 'Color',[0.7 0.7 0.7])
                        
                    end
                end
            end
            
            plot([0 1],[mean(group(grNo).first_pc_end(group(grNo).mouse_numbers)) mean(group(grNo).last_pc_start(group(grNo).mouse_numbers))],'-k','LineWidth', 3)
            CI = bootci(1000, @mean, group(grNo).first_pc_end(group(grNo).mouse_numbers));
            plot([0 0],CI,'-b','LineWidth',3)
            plot(0,mean(group(grNo).first_pc_end(group(grNo).mouse_numbers)),'ob','MarkerSize', 10,'MarkerFace','b')
            CI = bootci(1000, @mean, group(grNo).last_pc_start(group(grNo).mouse_numbers));
            plot([1 1],CI,'-r','LineWidth',3)
            plot(1,mean(group(grNo).last_pc_start(group(grNo).mouse_numbers)),'or','MarkerSize', 10,'MarkerFace','r')
            ylabel('Percent correct')
            ylim([30 110])
            title(['Percent correct, last ' num2str(trial_window) ' trials, pre and post'])
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            
            p_perCorr=ranksum(group(grNo).last_pc_start(group(grNo).mouse_numbers),group(grNo).first_pc_end(group(grNo).mouse_numbers));
            fprintf(1, '\np value for group No %d for ranksum test for percent correct= %d\n\n',grNo,p_perCorr);
            
            [h p_perCorrt]=ttest(group(grNo).last_pc_start(group(grNo).mouse_numbers),group(grNo).first_pc_end(group(grNo).mouse_numbers));
            fprintf(1, '\np value for group No %d for paired t test for percent correct= %d\n\n',grNo,p_perCorrt);
        end
        save([handles.drgb.outPathName '/' handles.drgb.outFileName],'group');
end






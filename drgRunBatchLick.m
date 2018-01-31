function drgRunBatchLick

%Ask user for the m file that contains information on what the user wants the analysis to be
%This file has all the information on what the user wants done, which files
%to process, what groups they fall into, etc
%
% An example of this file: drgbChoicesDanielPrelim
%
%

close all
clear all

trial_window=30;

lick_timecourseEv1=[];
lick_timecourseEv2=[];
inter_lick_intervals_Ev1=[];
inter_lick_intervals_Ev2=[];

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
    handles.time_end=5.2;
    handles.time_start=-2.2;
    
    if handles.lastTrialNo>=trial_window
        handles.trialNo=handles.lastTrialNo-trial_window;
        
        handles.evTypeNo=5;
        lick_timecourseEv1_this_file=[];
        inter_lick_intervalsEv1_this_file=[];
        [per_corr_per_trialEv1, lick_timecourseEv1_this_file,which_eventEv1,no_trialsEv1,inter_lick_intervalsEv1_this_file]=drgLickTimecourseThisEv(handles);
        lick_timecourseEv1=[lick_timecourseEv1; lick_timecourseEv1_this_file];
        inter_lick_intervals_Ev1=[inter_lick_intervals_Ev1 inter_lick_intervalsEv1_this_file];
        
        handles.evTypeNo=11;
        lick_timecourseEv2_this_file=[];
        inter_lick_intervalsEv2_this_file=[];
        [per_corr_per_trialEv2, lick_timecourseEv2_this_file,which_eventEv2,no_trialsEv2,inter_lick_intervalsEv2_this_file]=drgLickTimecourseThisEv(handles);
        lick_timecourseEv2=[lick_timecourseEv2; lick_timecourseEv2_this_file];
        inter_lick_intervals_Ev2=[inter_lick_intervals_Ev2 inter_lick_intervalsEv2_this_file];
        
    end
    
end

figure(1)
edges=[0:0.005:0.2];
h1=histogram(inter_lick_intervals_Ev1,edges)
h1.FaceColor='r';

hold on
h2=histogram(inter_lick_intervals_Ev2,edges)
h2.FaceColor='b';

plot([handles.smallest_inter_lick_interval handles.smallest_inter_lick_interval],[0 max([max(h1.Values) max(h2.Values)])],'r')
title('Histogram of the inter lick interval for S+ (red) and S- (blue)')
set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
                                    

%Timecourse for licks
try
    close 2
catch
end

hFig2 = figure(2);

no_lick_dt=floor(((handles.time_end-handles.time_pad)-(handles.time_start+handles.time_pad))/handles.dt_lick);
time=[1:no_lick_dt]*handles.dt_lick+handles.time_start+handles.time_pad;

Ev1_mean_licks=mean(lick_timecourseEv1,1)/handles.dt_lick;
Ev1_CI_licks = bootci(1000, @mean, lick_timecourseEv1)/handles.dt_lick;
mean_Ev1=[Ev1_mean_licks;Ev1_mean_licks];
Ev1_ci_l=Ev1_CI_licks(1,:)-mean_Ev1;

[hl1, hp1] = boundedline(time',Ev1_mean_licks', Ev1_ci_l', 'r');

hold on

Ev2_mean_licks=mean(lick_timecourseEv2,1)/handles.dt_lick;
Ev2_CI_licks = bootci(1000, @mean, lick_timecourseEv2)/handles.dt_lick;
mean_Ev2=[Ev2_mean_licks;Ev2_mean_licks];
Ev2_ci_l=Ev2_CI_licks(1,:)-mean_Ev2;

[hl2, hp2] = boundedline(time',Ev2_mean_licks', Ev2_ci_l', 'b');

title('Lick timecourse, S+ red, S- blue')
xlabel('Time (s)')
ylabel('Licks/s')
ylim([0 11.5])
set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)

pffft=1;
 
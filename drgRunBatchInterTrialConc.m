function drgRunBatchInterTrialConc

%This file calculates the estimated concentration based on the intertrial
%interval and the PID measurements of 20180618_sharpie_spmc_PID_180618_160444_out.mat
%Ask user for the m file that contains information on what the user wants the analysis to be
%This file has all the information on what the user wants done, which files
%to process, what groups they fall into, etc
%
% An example of this file: drgbChoicesDanielPrelim
%
%

close all
clear all

load('/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/pid/20180618_sharpie_spmc_PID_180618_160444_out.mat')
this_time=200000;
norm_fact=(f_before.a/2)*(1-exp(-this_time/f_before.b)+1-exp(-this_time/f_before.c));
t=0;
norm_fact_decay=(f_after.a/3)*(exp(-t/f_after.b)+exp(-t/f_after.c) +exp(-t/f_after.d) );

[choiceFileName,choiceBatchPathName] = uigetfile({'drgbChoices*.m'},'Select the .m file with all the choices for analysis');
addpath(choiceBatchPathName)
eval(['handles=' choiceFileName(1:end-2) ';'])
handles.choiceFileName=choiceFileName;
handles.choiceBatchPathName=choiceBatchPathName;

new_no_files=handles.drgbchoices.no_files;
choicePathName=handles.drgbchoices.PathName;
choiceFileName=handles.drgbchoices.FileName;

concIDs=[16 17 18 19 20 21 22];
fwd_concs=[10 3.2 1 0.32 0.1 0.032];
rev_concs=[0.032 0.1 0.32 1 3.2 10];

%For forward 16 is 10% and 22 is 0.032%
%For reversed 16 is 0.032% and 22 is 10%


%Very, very important!
handles.evTypeNo=2; %We do all OdorOn trials
LFPNo=1; %All trialss where LFP 1 can be read are processed

%If you want to skip files that have already been processed enter the number of the first file
% first_file=handles.drgb.first_file;

% if first_file==1
%     handles.drgb.lfpevpair_no=0;
%     handles.drgb.lfp_per_exp_no=0;
% else
%     load([handles.drgb.outPathName handles.drgb.outFileName])
%     handles.drgb=handles_drgb.drgb;
%     %The user may add new files
%     handles.drgbchoices.no_files=new_no_files;
%     handles.drgbchoices.PathName=choicePathName;
%     handles.drgbchoices.FileName=choiceFileName;
% end
%
% test_batch=handles.drgbchoices.test_batch;

%Do batch processing for each file
actual_concs_per_trial=[]
nominal_concs_per_trial=[];
inter_trial_intervals=[];
all_trial_no=0;
all_iti_no=0;

for filNum=1:handles.drgbchoices.no_files
    
    file_no=filNum
    
    %     if test_batch==1
    %         if handles.drgbchoices.group_no(filNum)==1
    %             handles.data_vs_simulate=5;
    %         else
    %             handles.data_vs_simulate=6;
    %         end
    %     end
    
    %read the jt_times files
    handles.jtfullName=[handles.drgbchoices.PathName,handles.drgbchoices.FileName{filNum}];
    handles.jtFileName=handles.drgbchoices.FileName{filNum};
    handles.jtPathName=handles.drgbchoices.PathName;
    
    
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
    
    %Now calculate the concentration for each trial
    %Enter trials
    firstTr=handles.trialNo;
    lastTr=handles.lastTrialNo;
    sessionNo=handles.sessionNo;
    
    no_trials=0;
    nominal_conc_per_trial_file=[];
    time_per_trial=[];
    
    for trNo=firstTr:lastTr
        
        
        evNo = drgFindEvNo(handles,trNo,sessionNo);
        if evNo~=-1
            
            
            excludeTrial=drgExcludeTrialLFP(handles.drg,handles.peakLFPNo,handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);
            
            if excludeTrial==0
                
                %Note: handles.peakLFPNo is the reference LFP
                [ref, trialNo, can_read1] = drgGetTrialLFPData(handles, LFPNo, evNo, handles.evTypeNo, handles.time_start+handles.time_pad, handles.time_end-handles.time_pad);
                
                if (can_read1==1)
                    
                    no_trials=no_trials+1;
                    time_per_trial(no_trials)=handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo);
                    
                    this_ii=0;
                    for ii=[16 17 18 19 20 21]
                        handles.evTypeNo=ii;
                        evNo = drgFindEvNo(handles,trNo,sessionNo);
                        if evNo~=-1
                            this_ii=ii;
                        end
                        
                    end
                    handles.evTypeNo=2;
                    
                  
                    this_jj=find(concIDs==this_ii);
                    if handles.drgbchoices.group_no(filNum)==1
                        %forward
                        this_nominal_conc=fwd_concs(this_jj);
                    else
                        %forward
                        this_nominal_conc=rev_concs(this_jj);
                    end
                    
                    nominal_conc_per_trial_file(no_trials)=this_nominal_conc;
                    
                    all_trial_no=all_trial_no+1;
                    nominal_concs_per_trial(all_trial_no)=this_nominal_conc;
                    
                    %Now calculate the concentration at peak
                    this_time=2.83;
                    nominal_conc_at_end=this_nominal_conc*(f_before.a/2)*(1-exp(-this_time/f_before.b)+1-exp(-this_time/f_before.c))/norm_fact;
                    if no_trials==1
                        actual_concs_per_trial(all_trial_no)=nominal_conc_at_end;
                        actual_conc_per_trial_file(no_trials)=nominal_conc_at_end;
                    else
                        delta_t=time_per_trial(no_trials)-time_per_trial(no_trials-1);
                        all_iti_no=all_iti_no+1;
                        inter_trial_intervals(all_iti_no)=delta_t;
                        this_decay_factor=(f_after.a/3)*(exp(-delta_t/f_after.b)+exp(-delta_t/f_after.c) +exp(-delta_t/f_after.d) )/norm_fact_decay;
                        residual_conc_last_trial=actual_conc_per_trial_file(no_trials-1)*this_decay_factor; 
                        actual_concs_per_trial(all_trial_no)=nominal_conc_at_end+residual_conc_last_trial;
                        actual_conc_per_trial_file(no_trials)=nominal_conc_at_end+residual_conc_last_trial;
                    end

                    
                end
            end
        end
        
        %end %if eventstamps...
    end %for evNo
    
end

try
    close 1
catch
end

hFig2 = figure(1);
edges=[0:5:120];
h1=histogram(inter_trial_intervals,edges)
title('Histogram of the inter trial intervals')

save([choicePathName 'odorant_spillover.mat'],'actual_concs_per_trial','nominal_concs_per_trial','inter_trial_intervals')
   
pffft=1;

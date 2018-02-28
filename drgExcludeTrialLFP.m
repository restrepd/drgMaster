<<<<<<< HEAD
function do_exclude=drgExcludeTrialLFP(drg,lfpno,event_time,sessionNo)
%function do_exclude=drsExcludeTrialLFP(dj,channel,event_time)
% this function outputs 1 if trial should be excluded

%channel=handles.lfp.lfpno;
do_exclude=0;

%Exclusions should only be done for LFP, not for sniff et al
if lfpno<=16
    
    %Find whether the trial was excluded
    
    trNo=find((event_time>=drg.session(sessionNo).start_times)&(event_time<=(drg.session(sessionNo).start_times+drg.session(sessionNo).sec_per_trial)));
    
    if  drg.session(sessionNo).trials_processed(trNo)==0
        do_exclude=1;
    end
    if drg.session(sessionNo).trials_x_channel_processed(lfpno,trNo)==0
        do_exclude=1;
    end
end

end





=======
function do_exclude=drgExcludeTrialLFP(drg,lfpno,event_time,sessionNo)
%function do_exclude=drsExcludeTrialLFP(dj,channel,event_time)
% this function outputs 1 if trial should be excluded

%channel=handles.lfp.lfpno;
do_exclude=0;

%Exclusions should only be done for LFP, not for sniff et al
if lfpno<=16
    
    %Find whether the trial was excluded
    
    trNo=find((event_time>=drg.session(sessionNo).start_times)&(event_time<=(drg.session(sessionNo).start_times+drg.session(sessionNo).sec_per_trial)));
    
    if  drg.session(sessionNo).trials_processed(trNo)==0
        do_exclude=1;
    end
    if drg.session(sessionNo).trials_x_channel_processed(lfpno,trNo)==0
        do_exclude=1;
    end
end

end





>>>>>>> 362cddba4db155729b2e4c347f0e0147096a5855

<<<<<<< HEAD
function do_exclude=drgExcludeTrial(drg,channel,event_time,sessionNo)
% this function outputs 1 if trial should be excluded

do_exclude=0;

%Find whether the trial was excluded
trNo=find((event_time>=drg.session(sessionNo).start_times)&(event_time<=(drg.session(sessionNo).start_times+drg.session(sessionNo).sec_per_trial)));

if  drg.session(sessionNo).trials_processed(trNo)==0
    do_exclude=1;
end
for lfpno=channel*4-3:channel*4
    if drg.session(sessionNo).trials_x_channel_processed(lfpno,trNo)==0
        do_exclude=1;
    end
end








=======
function do_exclude=drgExcludeTrial(drg,channel,event_time,sessionNo)
% this function outputs 1 if trial should be excluded

do_exclude=0;

%Find whether the trial was excluded
trNo=find((event_time>=drg.session(sessionNo).start_times)&(event_time<=(drg.session(sessionNo).start_times+drg.session(sessionNo).sec_per_trial)));

if  drg.session(sessionNo).trials_processed(trNo)==0
    do_exclude=1;
end
for lfpno=channel*4-3:channel*4
    if drg.session(sessionNo).trials_x_channel_processed(lfpno,trNo)==0
        do_exclude=1;
    end
end








>>>>>>> 362cddba4db155729b2e4c347f0e0147096a5855

function trNo = drgFindTrNo(handles,evNo,sessionNo,evTypeNo)
%
%   Finds the trial number corresponding to the evNo

if nargin==3
    evTypeNo=handles.evTypeNo;
end


if handles.drg.draq_p.sec_before_trigger==6
    %This is the code that was used for Justin's paper, it does not work when
    %longer trials are read out by drta
    [mintr, trNo]=min(abs((handles.drg.session(sessionNo).trial_start-handles.drg.session(sessionNo).events(evTypeNo).times(evNo))));
else
    [mintr, trNo]=min(abs(((handles.drg.session(sessionNo).trial_start+handles.drg.draq_p.sec_before_trigger+0.5)-handles.drg.session(sessionNo).events(evTypeNo).times(evNo))));
end




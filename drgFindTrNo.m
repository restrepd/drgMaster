function trNo = drgFindTrNo(handles,evNo,sessionNo,evTypeNo)
%
%   Finds the trial number corresponding to the evNo

if nargin==3
    evTypeNo=handles.evTypeNo;
end

[mintr, trNo]=min(abs((handles.drg.session(sessionNo).trial_start-handles.drg.session(sessionNo).events(evTypeNo).times(evNo))));




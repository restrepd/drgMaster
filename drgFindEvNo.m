function evNo = drgFindEvNo(handles,trNo,sessionNo,evTypeNo)
%
%   Finds the event number corresponding to the trNo

if nargin==3
    evTypeNo=handles.evTypeNo;
end

evNos=find ((handles.drg.session(sessionNo).events(evTypeNo).times>=handles.drg.session(sessionNo).trial_start(trNo))&...
    (handles.drg.session(sessionNo).events(evTypeNo).times<=(handles.drg.session(sessionNo).trial_start(trNo))...
    +handles.drg.session(sessionNo).sec_per_trial));

if isempty(evNos)
    evNo=-1;
else
    evNo=evNos(1);
end

end


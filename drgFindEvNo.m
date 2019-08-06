function evNo = drgFindEvNo(handles,trNo,sessionNo,evTypeNo)
%
%   Finds the event number corresponding to the trNo

if nargin==3
    evTypeNo=handles.evTypeNo;
end
%
% evNos=find ((handles.drg.session(sessionNo).events(evTypeNo).times>=handles.drg.session(sessionNo).trial_start(trNo))&...
%     (handles.drg.session(sessionNo).events(evTypeNo).times<=(handles.drg.session(sessionNo).trial_start(trNo))...
%     +handles.drg.session(sessionNo).sec_per_trial));
%
% if isempty(evNos)
%     evNo=-1;
% else
%     evNo=evNos(1);
% end

if handles.drg.draq_p.sec_before_trigger==6
    %This is the code that was used for Justin's paper, it does not work when
    %longer trials are read out by drta
    
    [this_min,evNo]=min(abs(handles.drg.session(sessionNo).events(evTypeNo).times-handles.drg.session(sessionNo).trial_start(trNo)));
    if this_min>handles.drg.session(sessionNo).draq_p.sec_before_trigger
        evNo=-1;
    end
    
    % %Find whether multiple trials match to this event and vet one out
    if evNo~=-1
        ii_trNo=find(abs(handles.drg.session(sessionNo).events(evTypeNo).times(evNo)-handles.drg.session(sessionNo).trial_start)<=handles.drg.session(sessionNo).draq_p.sec_before_trigger);
        if length(ii_trNo>1)
            delta_t=zeros(1,length(ii_trNo)-1);
            ii_delta_t=0;
            for ii=1:length(ii_trNo)
                if ii_trNo(ii)~=trNo
                    ii_delta_t=ii_delta_t+1;
                    delta_t(ii_delta_t)=abs(handles.drg.session(sessionNo).events(evTypeNo).times(evNo)-handles.drg.session(sessionNo).trial_start(ii_trNo(ii)));
                end
            end
            if sum(abs(handles.drg.session(sessionNo).events(evTypeNo).times(evNo)-handles.drg.session(sessionNo).trial_start(trNo))>...
                    delta_t)>0
                evNo=-1;
            end
        end
    end
    
else
    [this_min,evNo]=min(abs(handles.drg.session(sessionNo).events(evTypeNo).times-(handles.drg.session(sessionNo).trial_start(trNo)+handles.drg.draq_p.sec_before_trigger+0.5)));
    if this_min>5
        evNo=-1;
    end 
end

pffft=1;


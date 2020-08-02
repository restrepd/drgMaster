function [perCorr, encoding_trials, retrieval_trials,encoding_this_evTypeNo,retrieval_this_evTypeNo]=drgFindEncRetr(handles)
%This finds the encoding and retrieval trials for the current session
%This must be an spm session


sliding_window=20; %Trials for determination of behavioral performance
min_precent_high_beh=80; %Minimum percent correct for good behavior blocks
max_percent_low_beh=65;


sessionNo=handles.drg.unit(handles.unitNo).sessionNo;

switch handles.drg.drta_p.which_c_program
    case {2,10,14}
        %These are spm, etc
        odorOn=2;
        splus=5;
        hit=3;
        miss=7;
        CR=9;
        FA=13;
        sminus=11;
        if handles.drg.session(sessionNo).events(odorOn).noTimes-sliding_window+1<1
            sliding_window=handles.drg.session(sessionNo).events(odorOn).noTimes;
        end
        
    for ii=1:handles.drg.session(sessionNo).events(odorOn).noTimes-sliding_window+1
        first_time=handles.drg.session(sessionNo).events(odorOn).times(ii);
        last_time=handles.drg.session(sessionNo).events(odorOn).times(ii+sliding_window-1);
        noCRs=sum((handles.drg.session(sessionNo).events(CR).times>=first_time)&...
            (handles.drg.session(sessionNo).events(CR).times<=last_time));
        noHits=sum((handles.drg.session(sessionNo).events(hit).times>=first_time)&...
            (handles.drg.session(sessionNo).events(hit).times<=last_time));
        perCorr(ii+floor(sliding_window/2))=100*(noCRs+noHits)/sliding_window;

        if ii==1
            perCorr(ii:floor(sliding_window/2))=perCorr(ii+floor(sliding_window/2));
        end
        if ii==handles.drg.session(sessionNo).events(odorOn).noTimes-sliding_window+1
            perCorr(ii+floor(sliding_window/2)+1:handles.drg.session(sessionNo).events(odorOn).noTimes)=perCorr(ii+floor(sliding_window/2));
        end
    end
    
    
    encoding_trials=(perCorr<=max_percent_low_beh);
    retrieval_trials=(perCorr>=min_precent_high_beh);
    
    %Trials should not be encoding trials after the animal learns
    first_retrieval=find(retrieval_trials==1,1,'first');
    encoding_trials(first_retrieval:end)=0;
    
    %Now find the encoding and retreival event numbers for the current event
    %type
    encoding_ii=0;
    retrieval_ii=0;
    encoding_this_evTypeNo=[];
    retrieval_this_evTypeNo=[];
    for ii=1:handles.drg.session(sessionNo).events(odorOn).noTimes
        minii=[];
        [mindt minii]=min(abs(handles.drg.session(sessionNo).events(odorOn).times(ii)...
            -handles.drg.session(sessionNo).events(handles.evTypeNo).times));
        if mindt<handles.drg.session(sessionNo).draq_p.sec_before_trigger
            if encoding_trials(ii)==1
                encoding_ii=encoding_ii+1;
                encoding_this_evTypeNo(encoding_ii)=minii(1);
            end
            if retrieval_trials(ii)==1
                retrieval_ii=retrieval_ii+1;
                retrieval_this_evTypeNo(retrieval_ii)=minii(1);
            end
        end
        
    end
    
    otherwise
    %Not spm
    perCorr=100*ones(1,handles.drg.draq_d.noTrials);
    encoding_trials=ones(1,handles.drg.draq_d.noTrials);
    retrieval_trials=zeros(1,handles.drg.draq_d.noTrials);
    encoding_this_evTypeNo=1;
    retrieval_this_evTypeNo=1;
end







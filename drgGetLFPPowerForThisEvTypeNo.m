function [out_times,f,all_Power, all_Power_ref, all_Power_timecourse, this_trialNo, perCorr_pertr, which_event]=drgGetLFPPowerForThisEvTypeNo(handles)

%Generates a trial per trial phase histogram
odorOn=2;
sessionNo=handles.sessionNo;
Fs=handles.drg.session(sessionNo).draq_p.ActualRate;
freq=handles.burstLowF:(handles.burstHighF-handles.burstLowF)/100:handles.burstHighF;

window=round(handles.window*handles.drg.draq_p.ActualRate);
noverlap=round(handles.noverlap*handles.drg.draq_p.ActualRate);



%Enter trials
firstTr=handles.trialNo;
lastTr=handles.lastTrialNo;

no_trials=0;
all_Power=[];
all_Power_ref=[];
out_times=[];
f=[];
all_Power_timecourse=[];
this_trialNo=[];
perCorr_pertr=[];
which_event=[];
no_excluded=0;

% notch60HzFilt = designfilt('bandstopiir','FilterOrder',2, ...
%     'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
%     'DesignMethod','butter','SampleRate',floor(handles.drg.session(sessionNo).draq_p.ActualRate));

[perCorr, encoding_trials, retrieval_trials, encoding_this_evTypeNo,retrieval_this_evTypeNo]=drgFindEncRetr(handles);

for trNo=firstTr:lastTr
    
    if handles.save_drgb==0
        trialNo=trNo
    end
    
    evNo = drgFindEvNo(handles,trNo,sessionNo);
    
    if evNo~=-1
        
        
        excludeTrial=drgExcludeTrialLFP(handles.drg,handles.peakLFPNo,handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);
        
        
        if excludeTrial==0
            
            %First find the time range for the spectrogram
            if handles.subtractRef==1
                if handles.time_start+handles.time_pad<handles.startRef+handles.time_pad
                    min_t=handles.time_start+handles.time_pad;
                else
                    min_t=handles.startRef+handles.time_pad;
                end
                
                if handles.time_end-handles.time_pad>handles.endRef-handles.time_pad
                    max_t=handles.time_end-handles.time_pad;
                else
                    max_t=handles.endRef-handles.time_pad;
                end
            else
                min_t=handles.time_start+handles.time_pad;
                max_t=handles.time_end-handles.time_pad;
            end
            
            
            [LFP, trialNo, can_read] = drgGetTrialLFPData(handles, handles.peakLFPNo, evNo, handles.evTypeNo, min_t, max_t);
            
            if (can_read==1)
                
                
                %Note: I tried hamming, hann and flattopwin widows, the
                %results are qualitatively different, but they look the
                %same
                
                
                [S,f,t,P]=spectrogram(detrend(double(LFP)),window,noverlap,freq,handles.drg.session(handles.sessionNo).draq_p.ActualRate);
                
                no_trials=no_trials+1;
                this_trialNo(no_trials)=trNo;
                
                times=t+min_t;
                out_times=times((times>=handles.time_start+handles.time_pad)&(times<=handles.time_end-handles.time_pad));
                all_Power_timecourse(no_trials,1:length(f),1:length(out_times))=P(:,(times>=handles.time_start+handles.time_pad)&(times<=handles.time_end-handles.time_pad));
                all_Power(no_trials,1:length(f))=mean(P((times>=handles.time_start+handles.time_pad)&(times<=handles.time_end-handles.time_pad)),2);
                if handles.subtractRef==1
                    P_ref=[];
                    P_ref=P(:,(times>=handles.startRef+handles.time_pad)&(times<=handles.endRef-handles.time_pad));
                    all_Power_ref(no_trials,1:length(f))=mean(P_ref,2);
                end
                switch handles.drg.drta_p.which_c_program
                    case {2,10}
                        perCorr_pertr(no_trials)=perCorr(drgFindEvNo(handles,trialNo,sessionNo,odorOn));
                    otherwise
                        perCorr_pertr(no_trials)=100;
                end
                if handles.displayData==0
                    for evTypeNo=1:length(handles.drgbchoices.evTypeNos)
                        switch handles.evTypeNo
                            case 1
                                %tstart is the reference event
                                if handles.drgbchoices.evTypeNos(evTypeNo)==1
                                    %This is tstart
                                    if sum(handles.drg.session(1).events(handles.drgbchoices.evTypeNos(evTypeNo)).times==handles.drg.session(1).events(handles.drgbchoices.referenceEvent).times(evNo))>0
                                        which_event(evTypeNo,no_trials)=1;
                                    else
                                        which_event(evTypeNo,no_trials)=0;
                                    end
                                else
                                    %These are not tstart, and the time
                                    %should be compared at OdorOn
                                    %This is tstart
                                    if sum(handles.drg.session(1).events(handles.drgbchoices.evTypeNos(evTypeNo)).times==handles.drg.session(1).events(2).times(evNo))>0
                                        which_event(evTypeNo,no_trials)=1;
                                    else
                                        which_event(evTypeNo,no_trials)=0;
                                    end
                                end
                            otherwise
                                %OdorOn is the reference event
                                if sum(handles.drg.session(1).events(handles.drgbchoices.evTypeNos(evTypeNo)).times==handles.drg.session(1).events(handles.drgbchoices.referenceEvent).times(evNo))>0
                                    which_event(evTypeNo,no_trials)=1;
                                else
                                    which_event(evTypeNo,no_trials)=0;
                                end
                        end
                    end
                end
                
            end
        end
    else
        no_excluded=no_excluded+1;
    end %for evNo
    
end



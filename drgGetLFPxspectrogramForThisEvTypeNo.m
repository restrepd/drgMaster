function [out_times,freq, all_xspec_ref,all_xspec_timecourse, this_trialNo, perCorr_pertr, which_event]=drgGetLFPxspectrogramForThisEvTypeNo(handles)

%Generates a cross-spectrogram

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
out_times=[];
all_xspec_timecourse=[];
this_trialNo=[];
perCorr_pertr=[];
which_event=[];
no_excluded=0;

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
            
                
            [LFP1, trialNo1, can_read1] = drgGetTrialLFPData(handles, handles.peakLFPNo, evNo, handles.evTypeNo, min_t, max_t);
            [LFP2, trialNo2, can_read2] = drgGetTrialLFPData(handles, handles.burstLFPNo, evNo, handles.evTypeNo, min_t, max_t);
            
            if (can_read1==1)&(can_read2==1)
                
                
                
                no_time_points=floor(length(LFP1)/(window-noverlap))-(window/(window-noverlap))+1;
                
                %Compute the crosspectrogram
                [xspec,f,txspec] = xspectrogram(LFP1,LFP2,hamming(window),noverlap,freq,handles.drg.session(handles.sessionNo).draq_p.ActualRate,'mimo');
                %note: mesh(t,f,20*log10(s))
                
                no_trials=no_trials+1;
                this_trialNo(no_trials)=trNo;
                
                out_times=txspec+min_t;
                
                these_xspecs=zeros(1,length(freq),length(out_times));
                these_xspecs(1,:,:)=xspec;
                all_xspec_timecourse(no_trials,1:length(freq),1:length(out_times))=xspec;
                
                 if handles.subtractRef==1
                    xspec_ref=[];
                    xspec_ref=xspec(:,(out_times>=handles.startRef+handles.time_pad)&(out_times<=handles.endRef-handles.time_pad));
                    all_xspec_ref(no_trials,1:length(f))=mean(xspec_ref,2);
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
    pffft=1;
end


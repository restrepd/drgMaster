function [out_times,freq,all_Cxy_timecourse, trial_numbers, perCorr_pertr, which_event]=drgGetLFPCoherenceForThisEvTypeNo(handles)

%Generates a coherence timecourse
%Borjigin et al. www.pnas.org/cgi/content/short/1308285110
%Coherence between EEG channels is measured by amplitude squared coherence 
%Cxy(f) (mscohere.m in MATLAB signal toolbox; MathWorks, Inc.), which is a 
%coherence estimate of the input signals x and y using Welch?s averaged, 
%modified periodogram method. The magnitude squared coherence estimate 
%Cxy(f) is a function of frequency with values between 0 and 1 that 
%indicates how well x corresponds to y at each frequency.

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
all_Cxy_timecourse=[];
trial_numbers=[];
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
            
            min_t=handles.time_start+handles.time_pad;
            max_t=handles.time_end-handles.time_pad;

                
            [LFP1, trialNo, can_read1] = drgGetTrialLFPData(handles, handles.peakLFPNo, evNo, handles.evTypeNo, min_t, max_t);
            [LFP2, trialNo, can_read2] = drgGetTrialLFPData(handles, handles.burstLFPNo, evNo, handles.evTypeNo, min_t, max_t);
            
            if (can_read1==1)&(can_read2==1)
                
                
                %Get coherence spectrogram
                %Estimate Cxy
                current_ii=0;
                %                 no_Cxys=0;
                no_time_points=floor(length(LFP1)/(window-noverlap))-(window/(window-noverlap));
                if no_time_points<=0
                    no_time_points=1;
                end
                Cxy=zeros(length(freq),no_time_points);
                tCxy=zeros(no_time_points,1);
                %while current_ii+window<length(LFP1)
                for no_Cxys=1:no_time_points
                    current_ii=(no_Cxys-1)*(window-noverlap);
                    tCxy(no_Cxys,1)=(current_ii+window/2)/handles.drg.session(handles.sessionNo).draq_p.ActualRate;
                    [thisCxy,fCxy] = mscohere(LFP1(current_ii+1:current_ii+window),LFP2(current_ii+1:current_ii+window),hann(window/2),noverlap/2,freq,handles.drg.session(handles.sessionNo).draq_p.ActualRate,'mimo');
                    Cxy(1:length(freq),no_Cxys)=thisCxy';
                    %                     current_ii=current_ii+window-noverlap;
                end
                
                %                 mesh(tCxy,fCxy,Cxy)
                %                 view(2)
                %                 axis tight
                %
                %                 %Compute the crosspectrogram
                %                 [s,f,t] = xspectrogram(LFP1,LFP2,hamming(window),noverlap,freq,handles.drg.session(handles.sessionNo).draq_p.ActualRate,'mimo');
                %                 %note: mesh(t,f,20*log10(s))
                %
                no_trials=no_trials+1;
                trial_numbers(no_trials)=trNo;
                
                out_times=tCxy+min_t;
                
                these_Cxys=zeros(1,length(freq),length(out_times));
                these_Cxys(1,:,:)=Cxy;
                all_Cxy_timecourse(no_trials,1:length(freq),1:length(out_times))=Cxy;
                
                
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


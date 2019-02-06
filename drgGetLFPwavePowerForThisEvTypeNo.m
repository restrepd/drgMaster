function [out_t,f,all_Power, all_Power_ref, all_Power_timecourse, this_trialNo, perCorr_pertr, which_event]=drgGetLFPwavePowerForThisEvTypeNo(handles)

%Generates a trial per trial phase histogram
odorOn=2;
sessionNo=handles.sessionNo;
Fs=handles.drg.session(sessionNo).draq_p.ActualRate;
dec_n=fix(handles.drg.session(sessionNo).draq_p.ActualRate/1000);
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

%Setup the wavelet scales
%   scales = helperCWTTimeFreqVector(minfreq,maxfreq,f0,dt,NumVoices)
%   f0 - center frequency of the wavelet in cycles/unit time
%   dt - sampling interval
%   NumVoices - number of voices per octave

NumVoices=5;
minfreq=handles.burstLowF;
maxfreq=handles.burstHighF;
dt=1/Fs;
f0=5/(2*pi);

a0 = 2^(1/NumVoices);
minscale = f0/(maxfreq*dt);
maxscale = f0/(minfreq*dt);
minscale = floor(NumVoices*log2(minscale));
maxscale = ceil(NumVoices*log2(maxscale));
scales = a0.^(minscale:maxscale).*dt;

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
            
            %First find the time range for the spectrogram, calculate using
            %handles.time_pad to exclude artifacts at the end of
            if handles.subtractRef==1
                if handles.time_start<handles.startRef
                    min_t=handles.time_start;
                else
                    min_t=handles.startRef;
                end
                
                if handles.time_end>handles.endRef
                    max_t=handles.time_end;
                else
                    max_t=handles.endRef;
                end
            else
                min_t=handles.time_start;
                max_t=handles.time_end;
            end
            

                
            [LFP, trialNo, can_read] = drgGetTrialLFPData(handles, handles.peakLFPNo, evNo, handles.evTypeNo, min_t, max_t);
%             load('/Users/restrepd/Documents/Projects/Ethan/MEA to dg/thisLFP.mat')
%             LFP(1,1:119281)=thisLFP(1,:);
            
            if (can_read==1)
                
   
                %Now do the wavelet transform
                decLFP=decimate(LFP,dec_n);
                decFs=Fs/dec_n;
               
                cwtLFP = cwtft({detrend(double(decLFP)),1/decFs},'wavelet','morl','scales',scales);
                Prev=abs(cwtLFP.cfs).^2;
                P=Prev(end:-1:1,:);
                DT=1/decFs;
                t = 0:DT:(numel(decLFP)*DT)-DT;
                frev=cwtLFP.frequencies;
                f=frev(end:-1:1);
                
                no_trials=no_trials+1;
                this_trialNo(no_trials)=trNo;
                
                
                all_times=t+min_t;
                out_times=all_times(:,(all_times>=handles.time_start+handles.time_pad)&(all_times<=handles.time_end-handles.time_pad));
                Pout(:,:)=P(:,(all_times>=handles.time_start+handles.time_pad)&(all_times<=handles.time_end-handles.time_pad));
                all_Power_timecourse(no_trials,1:length(f),1:length(out_times))=Pout(:,:);
                all_Power(no_trials,1:length(f))=mean(Pout,2);
                if handles.subtractRef==1
                    P_ref=[];
                    P_ref=P(:,(all_times>=handles.startRef+handles.time_pad)&(all_times<=handles.endRef-handles.time_pad));
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
        end %for evNo
        
    end
    
end

%Sometimes the number of points differs by one point between trials.
%Generate a time vector with the correct number of points
if no_trials>0
    sz_apt=size(all_Power_timecourse);
    t = 0:DT:(sz_apt(3)*DT)-DT;
    out_t=t+out_times(1);
else
    out_t=[];
end


function [out_times,freq,all_Cxy_timecourse, trial_numbers, perCorr_pertr, which_event]=drgGetLFPimgCoherenceForThisEvTypeNo(handles)

%Generates imaginary coherence timecourse
%Nolte et al Clinical Neuropsych 115:2292-2307, 2004
%https://stackoverflow.com/questions/41923700/imaginary-part-of-coherence-matlab
% Assuming the processes producing x and y are ergodic, the expectations can be estimated by computing the average over many blocks of data. With that in mind, an implementation of the coherency as described in your definition could look like:
%
% function [ result ] = coherency( x,y,N )
%   % divide data in N equal length blocks for averaging later on
%   L  = floor(length(x)/N);
%   xt = reshape(x(1:L*N), L, N);
%   yt = reshape(y(1:L*N), L, N);
%
%   % transform to frequency domain
%   Xf = fft(xt,L,1);
%   Yf = fft(yt,L,1);
%
%   % estimate expectations by taking the average over N blocks
%   xy = sum(Xf .* conj(Yf), 2)/N;
%   xx = sum(Xf .* conj(Xf), 2)/N;
%   yy = sum(Yf .* conj(Yf), 2)/N;
%
%   % combine terms to get final result
%   result=xy./sqrt(xx.*yy);
% end
% If you only want the imaginary part, then it's a simple matter of computing imag(coherency(x,y,N)).
 
odorOn=2;
sessionNo=handles.sessionNo;
Fs=handles.drg.session(sessionNo).draq_p.ActualRate;


window=round(handles.window*handles.drg.draq_p.ActualRate);
noverlap=round(handles.noverlap*handles.drg.draq_p.ActualRate);

%Resample four times within the window
N=4;
L  = floor(window/N);
f=Fs*(0:(L/2))/L;
freq=f((f>=handles.burstLowF)&(f<=handles.burstHighF));

 
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
            max_t=handles.time_end; %Note: I leave the pad at the end to prevent errors
            
            
            [LFP1, trialNo, can_read1] = drgGetTrialLFPData(handles, handles.peakLFPNo, evNo, handles.evTypeNo, min_t, max_t);
            [LFP2, trialNo, can_read2] = drgGetTrialLFPData(handles, handles.burstLFPNo, evNo, handles.evTypeNo, min_t, max_t);
            
            %Test the algorithm, in this example LFP1 has 10 Hz with an onset 20 msec
            %before LFP2
%             time=(1:length(LFP1))*(1/Fs);
%             this_f=10;
%             LFP1=4000*sin(2*pi()*this_f*time)+1000*randn(1,length(LFP1));
%             delay=0.02;
%             LFP2=4000*sin(2*pi()*this_f*(time+delay))+1000*randn(1,length(LFP1));
%             
%             figure(5)
%             plot(time(1:20000),LFP1(1:20000),'-b')
%             hold on
%             plot(time(1:20000),LFP2(1:20000),'-r')
%             legend('LFP1','LFP2')
%             title('10 Hz with LFP2 20 msec before LFP1')
            
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
                    %                     [thisCxy,fCxy] = mscohere(LFP1(current_ii+1:current_ii+window),LFP2(current_ii+1:current_ii+window),hann(window/2),noverlap/2,freq,handles.drg.session(handles.sessionNo).draq_p.ActualRate,'mimo');
                    
                    x=LFP1(current_ii+1:current_ii+window);
                    y=LFP2(current_ii+1:current_ii+window);
                    
                    % divide data in N equal length blocks for averaging later on
                    xt = reshape(x(1:L*N), L, N);
                    yt = reshape(y(1:L*N), L, N);
                    
                    % transform to frequency domain
                    Xf = fft(xt,L,1);
                    Yf = fft(yt,L,1);
                    
                    % estimate cross spectrum and auto spectra
                    xy = sum(Xf .* conj(Yf), 2)/N;
                    xx = sum(Xf .* conj(Xf), 2)/N;
                    yy = sum(Yf .* conj(Yf), 2)/N;
                    
                    % combine terms to get final result
                    this_imagCxy=imag(xy./sqrt(xx.*yy));
                    this_imagCxy=this_imagCxy((f>=handles.burstLowF)&(f<=handles.burstHighF));
                    
                    Cxy(1:length(freq),no_Cxys)=this_imagCxy';
                    
                    
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


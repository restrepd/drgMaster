function drgLFPspect(handles)

%Generates LFP spectrogram for all trials
sessionNo=handles.sessionNo;
Fs=handles.drg.session(sessionNo).draq_p.ActualRate;
freq=handles.burstLowF:ceil((handles.burstHighF-handles.burstLowF)/50):handles.burstHighF;
highF1=handles.burstLowF;
highF2=handles.burstHighF;
pad_time=handles.time_pad;
n_phase_bins=handles.n_phase_bins;

window=round(handles.window*handles.drg.draq_p.ActualRate); 
noverlap=round(handles.noverlap*handles.drg.draq_p.ActualRate); 




%Now get the phase of gamma bursts witin theta
no_trials=0;
no_ref_trials=0;
all_Power=[];
all_Power_ref=[];

firstTr=handles.trialNo;
lastTr=handles.lastTrialNo;


for trNo=firstTr:lastTr
    
    if handles.save_drgb==0
        trialNo=trNo
    end
    
    evNo = drgFindEvNo(handles,trNo,sessionNo);
    
    if evNo~=-1
        
        excludeTrial=drgExcludeTrialLFP(handles.drg,handles.peakLFPNo,handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);
        
        if excludeTrial==0
            
            
            [LFPhigh, trialNo, can_read] = drgGetTrialLFPData(handles, handles.peakLFPNo, evNo, handles.evTypeNo, handles.time_start+handles.time_pad, handles.time_end-handles.time_pad);
            if handles.subtractRef==1
                [LFPhigh_ref, trialNo, can_read_ref] = drgGetTrialLFPData(handles, handles.peakLFPNo, evNo, handles.evTypeNo, handles.startRef+handles.time_pad, handles.endRef-handles.time_pad);
                can_read=can_read&can_read_ref;
            end
            
            if (can_read==1)
                [S,f,t,P]=spectrogram(detrend(double(LFPhigh)),window,noverlap,freq,handles.drg.session(handles.sessionNo).draq_p.ActualRate);
                
                if handles.subtractRef==1
                    %Get the reference LFP
                    [S,f,t,P_ref]=spectrogram(detrend(double(LFPhigh_ref)),window,noverlap,freq,handles.drg.session(handles.sessionNo).draq_p.ActualRate);
                end
                
                %Sometimes the power is zero
                not_inf=1;
                if sum(isinf(mean(10*log10(P),2)))>0
                    not_inf=0;
                end
                
                
                if (handles.subtractRef==1)
                    if sum(isinf(mean(10*log10(P_ref),2)))>0
                        not_inf=0;
                    end
                end   

                if not_inf==1
                    no_trials=no_trials+1;
                    log_all_Power(no_trials,:)=mean(10*log10(P),2);
                    if (handles.subtractRef==1)
                        log_all_Power_ref(no_trials,:)=mean(10*log10(P_ref),2);
                    end
                end
                
            end
            
        end
        %end
        %end %if eventstamps...
    end %for evNo
end

% %Now plot the encoding log10P
try
    close 1
catch
end

hFig1 = figure(1);
set(hFig1, 'units','normalized','position',[.25 .25 .23 .65])

if handles.subtractRef==0
    shadedErrorBar(freq,mean(log_all_Power,1),std(log_all_Power,0,1)/sqrt(no_trials),'-b')
    hold on
    %errorbar(freq,mean(log10(all_Power),1),std(log10(all_Power),0,1)/sqrt(no_trials),'-b')
    title('Power vs frequency for all trials')
else
    shadedErrorBar(freq,mean(log_all_Power-log_all_Power_ref,1),std(log_all_Power-log_all_Power_ref,0,1)/sqrt(no_trials),'-b')
    hold on
    %errorbar(freq,mean(log10(all_Power)-log10(all_Power_ref),1),std(log10(all_Power)-log10(all_Power_ref),0,1)/sqrt(no_trials),'-b')
    title('Power vs frequency for all trials (reference subtracted)')
end
ylabel('Power (dB)')
xlabel('Frequency (Hz)')

if handles.subtractRef==0
    if handles.autoscale==1
        maxlp=prctile(mean(log_all_Power,1),99);
        minlp=prctile(mean(log_all_Power,1),1);
    else
        maxlp=handles.maxLogP;
        minlp=handles.minLogP;
    end
else
    
    if handles.autoscale==1
        maxlp=prctile(mean(log_all_Power,1)-mean(log_all_Power_ref,1),99);
        minlp=prctile(mean(log_all_Power,1)-mean(log_all_Power_ref,1),1);
    else
        maxlp=handles.maxLogP;
        minlp=handles.minLogP;
    end
    
end

%Note: Diego added this on purpose to limit the range to 10 dB
%This results in emphasizing changes in the top 10 dB
if maxlp-minlp>10
    minlp=maxlp-10;
end

try
    close 2
catch
end
 
hFig2 = figure(2);
set(hFig2, 'units','normalized','position',[.01 .1 .23 .8])

trials=1:no_trials;
if handles.subtractRef==0
    drg_pcolor(repmat(freq,length(trials),1),repmat(trials',1,length(freq)),log_all_Power)
else
    drg_pcolor(repmat(freq,length(trials),1),repmat(trials',1,length(freq)),log_all_Power-log_all_Power_ref)   
end
colormap jet
shading flat
caxis([minlp maxlp]);
xlabel('Frequency (Hz)')
ylabel('Trial No');
title(['log10(Power) per trial ' handles.drg.session(1).draq_d.eventlabels{handles.evTypeNo}])

try
    close 3
catch
end

hFig3 = figure(3);
set(hFig3, 'units','normalized','position',[.49 .1 .05 .3])

prain=[minlp:(maxlp-minlp)/99:maxlp];
drg_pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
colormap jet
shading interp
ax=gca;
set(ax,'XTickLabel','')

pffft=1


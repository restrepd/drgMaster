function handles=drgPLVTimecourse(handles)

if ~isfield(handles,'randpermLFP')
    handles.randpermLFP=0;
end

%Generates a trial per trial phase histogram
odorOn=2;
 
sessionNo=handles.sessionNo;
Fs=handles.drg.session(sessionNo).draq_p.ActualRate;
lowF1=handles.peakLowF;
lowF2=handles.peakHighF;
highF1=handles.burstLowF;
highF2=handles.burstHighF;
pad_time=handles.time_pad;
n_phase_bins=handles.n_phase_bins;
decimation_factor=40;

%Empty vectors
handles.drgb.LFPcorr=[];

%Enter trials
firstTr=handles.trialNo;
lastTr=handles.lastTrialNo;

[perCorr, encoding_trials, retrieval_trials, encoding_this_evTypeNo,retrieval_this_evTypeNo]=drgFindEncRetr(handles);

no_trials=0;
trials_attempted=0;

vetted_trNos=[];
vetted_evNos=[];
ii_vet=0;
these_LFPs=[];
for trNo=firstTr:lastTr
    evNo = drgFindEvNo(handles,trNo,sessionNo);
    if evNo~=-1
        excludeTrial=drgExcludeTrialLFP(handles.drg,handles.peakLFPNo,handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);
        
        if excludeTrial==0
            [LFP1, trialNo, can_read1] = drgGetTrialLFPData(handles, handles.peakLFPNo, evNo, handles.evTypeNo, handles.time_start, handles.time_end);
            
            [LFP2, trialNo, can_read2] = drgGetTrialLFPData(handles, handles.burstLFPNo, evNo, handles.evTypeNo, handles.time_start, handles.time_end);
            num_timepoints=floor(length(LFP1)/decimation_factor);
           
            
            if (can_read1==1)&(can_read2==1)
                ii_vet=ii_vet+1;
                handles.drgb.PLV.no_trials=ii_vet;
                handles.drgb.PLV.perCorr(ii_vet)=perCorr(drgFindEvNo(handles,trNo,sessionNo,odorOn));
                handles.drgb.PLV.trial(ii_vet).trialNo=trNo;
                if handles.save_drgb==1
                    for evTypeNo=1:length(handles.drgbchoices.evTypeNos)
                        if sum(handles.drg.session(1).events(handles.drgbchoices.evTypeNos(evTypeNo)).times==handles.drg.session(1).events(handles.drgbchoices.referenceEvent).times(evNo))>0
                            handles.drgb.PLV.which_event(evTypeNo,ii_vet)=1;
                        else
                            handles.drgb.PLV.which_event(evTypeNo,ii_vet)=0;
                        end
                    end
                end
                
                vetted_trNos(ii_vet)=trNo;
                vetted_evNos(ii_vet)=evNo;
                LFP1=decimate(LFP1,decimation_factor);
                these_LFPs(1,1:length(LFP1),ii_vet)=LFP1;
                LFP2=decimate(LFP2,decimation_factor);
                these_LFPs(2,1:length(LFP2),ii_vet)=LFP2;
                if handles.displayData==1
                    fprintf(1, ['Trial No %d, event No %d\n'], trNo,evNo);
                end 
            else
                if handles.displayData==1
                    fprintf(1, ['Trial No %d, event No %d, LFP could not be read\n'], trNo,evNo);
                end
            end
        else
            if handles.displayData==1
                fprintf(1, ['Trial No %d, event No %d, trial excluded\n'], trNo,evNo);
            end
        end
    else
        if handles.displayData==1
            fprintf(1, ['Trial No %d, event No %d\n'], trNo,evNo);
        end
    end
    %end
    %end %if eventstamps...
end %for evNo


%Calculate PLV
filtSpec.range = [lowF1 lowF2];
filtSpec.order=50;
[plv,delta_phase]=pn_eegPLV(these_LFPs, Fs/decimation_factor, filtSpec);

%Discard the ends
ii_start=floor(handles.time_pad*Fs/decimation_factor)+1;
ii_end=ii_start+floor((handles.time_end-handles.time_start-2*handles.time_pad)*Fs/decimation_factor)-1;
plv = plv(ii_start:ii_end,:,:);
delta_phase=delta_phase(ii_start:ii_end,:,:);
handles.drgb.PLV.plv=plv;
handles.drgb.PLV.delta_phase=delta_phase;

time=handles.time_start+handles.time_pad:1/floor(Fs/decimation_factor):handles.time_end-handles.time_pad;
time=time(1:length(plv));
handles.drgb.PLV.time=time;

%Get shuffled PLV
% if handles.randpermLFP==1
perm_ii_vet = randperm(length(vetted_trNos));
shuffled_LFPs=[];
shuffled_LFPs=these_LFPs;
for ii_vet=1:length(vetted_trNos)
    
    trNo_sh=vetted_trNos(perm_ii_vet(ii_vet));
    evNo_sh = drgFindEvNo(handles,trNo_sh,sessionNo);
    
    [LFP2, trialNo, can_read2] = drgGetTrialLFPData(handles, handles.burstLFPNo, evNo_sh, handles.evTypeNo, handles.time_start, handles.time_end);
    LFP2=decimate(LFP2,decimation_factor);
    shuffled_LFPs(2,:,ii_vet)=LFP2;
    if handles.displayData==1
        fprintf(1, ['Shuffled trial No %d, event No %d\n'], trNo_sh,evNo_sh);
    end
end %for evNo


%Calculate shuffled PLV
[plv_sh,delta_phase_sh]=pn_eegPLV(shuffled_LFPs, Fs/decimation_factor, filtSpec);

%Discard the ends
ii_start=floor(handles.time_pad*Fs/decimation_factor)+1;
ii_end=ii_start+floor((handles.time_end-handles.time_start-2*handles.time_pad)*Fs/decimation_factor)-1;
plv_sh = plv_sh(ii_start:ii_end,:,:);
handles.drgb.PLV.plv_sh=plv_sh;
delta_phase_sh = delta_phase_sh(ii_start:ii_end,:,:);
handles.drgb.PLV.delta_phase_sh=delta_phase_sh;
% end


%Calculate the weigted mean t_lag

if handles.displayData==1

    
    
    %Plot PLV timecourse

    try
        close 1
    catch
    end
    
    
    hFig1 = figure(1);
    set(hFig1, 'units','normalized','position',[.05 .33 .45 .25])
    
    this_plv=plv(:,1,2);
    this_plv=this_plv';
    plot(time,this_plv,'-b')
    
    hold on
    
    pls=prctile(plv_sh(:,1,2),95);
    plot([time(1) time(end)],[pls pls],'-r')
    
    xlabel('Time (sec)')
    ylabel('PLV')
    title('PLV timecourse')
    
    %Plot PLV trajectory in polar coordinates
    try
        close 2
    catch
    end
    
    hFig2 = figure(2);
    set(hFig2, 'units','normalized','position',[.05 .33 .45 .25])
    
    this_delta_phase=delta_phase(:,1,2);
    this_delta_phase=this_delta_phase';
    
    p = polarplot(this_delta_phase,this_plv);
    title('PLV timecourse polar plot')
    
   %Plot phase timecourse
    try
        close 3
    catch
    end
    
    hFig3 = figure(3);
    set(hFig3, 'units','normalized','position',[.05 .33 .45 .25])
    
    plot(time,this_delta_phase,'-b')
     xlabel('Time (sec)')
    ylabel('Delta phase (rad)')
    title('Delta phase timecourse')
end
     
pfffft=1;






function handles=drgThetaAmpPhaseTrialRangeConc(handles)

%Generates a per trial phase histogram for the concentration experimens of
%Justin
%Note: handles.evTypeNo is 16 for HiOd1 and 21 for LowOd6 forward

odorOn=2;
splus=5;
first_conc_evNo=16;
last_conc_evNo=21;
sessionNo=handles.sessionNo;
Fs=handles.drg.session(sessionNo).draq_p.ActualRate;
lowF1=handles.peakLowF;
lowF2=handles.peakHighF;
highF1=handles.burstLowF;
highF2=handles.burstHighF;
pad_time=handles.time_pad;
n_phase_bins=handles.n_phase_bins;

%Empty vectors
handles.drgb.PAC.no_trials=0;
handles.drgb.PAC.meanVectorLength=[];
handles.drgb.PAC.meanVectorAngle=[];
handles.drgb.PAC.peakAngle=[];
handles.drgb.PAC.mod_indx=[];
handles.drgb.PAC.this_trialNo=[];
handles.drgb.PAC.all_phase_histo=[];
handles.drgb.PAC.all_theta_wave=[];
handles.drgb.PAC.perCorr=[];
handles.drgb.PAC.which_event=[];

%Enter trials
firstTr=handles.trialNo;
lastTr=handles.lastTrialNo;


[perCorr, encoding_trials, retrieval_trials, encoding_this_evTypeNo,retrieval_this_evTypeNo]=drgFindEncRetr(handles);
%handles=drgPercentLickPerTrial(handles);

% if handles.displayData==1
%     try
%         close 1
%     catch
%     end
%
%     hFig1 = figure(1);
%     set(hFig1, 'units','normalized','position',[.55 .27 .35 .15])
% end

 
%Now get the phase of gamma bursts witin theta
no_encoding_trials=0;
no_retrieval_trials=0;
no_trials=0;
enc_phase_histo=[];
retr_phase_histo=[];
all_phase_histo=[];
all_out_times=[];
all_theta_wave=[];
MI_enc=[];
MI_retr=[];
which_event=[];
all_out_time_PAChisto=[];
spm=[];
trials_attempted=0;
conc_evTyNoPerTr=[];
currentevTypeNo=handles.evTypeNo;

for conc_evTyNo=last_conc_evNo:-1:first_conc_evNo
    
    handles.evTypeNo=conc_evTyNo;
    
    for trNo=firstTr:lastTr
        
        evNo = drgFindEvNo(handles,trNo,sessionNo);
        
        if evNo~=-1
            excludeTrial=drgExcludeTrialLFP(handles.drg,handles.peakLFPNo,handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);
            
            if excludeTrial==0
                
                trials_attempted=trials_attempted+1;
                [LFPlow, trialNo, can_read1] = drgGetTrialLFPData(handles, handles.peakLFPNo, evNo, handles.evTypeNo, handles.time_start, handles.time_end);
                %This excludes flat lick recordings
                if handles.peakLFPNo==19
                    if sum(LFPlow)==0
                        can_read1=0;
                    end
                end
                [LFPhigh, trialNo, can_read2] = drgGetTrialLFPData(handles, handles.burstLFPNo, evNo, handles.evTypeNo, handles.time_start, handles.time_end);
                
                if (can_read1==1)&(can_read2==1)
                    
                    no_trials=no_trials+1;
                    
                    if handles.displayData==1
                        fprintf(1, ['Trial No %d, event No %d, processed trial no %d\n'], trNo,evNo,no_trials);
                    end
                    
                    %If this is being used for batch processing find out whether
                    %this is an S+
                    eventNo = drgFindEvNo(handles,trNo,sessionNo,splus);
                    if eventNo~=-1
                        spm(no_trials)=1;
                    else
                        spm(no_trials)=0;
                    end
                    
                    time(no_trials)=handles.drg.session(sessionNo).trial_start(trNo);
                    handles.drgb.PAC.this_trialNo(no_trials)=trNo;
                    if trNo==1
                        handles.drgb.PAC.delta_t_trial(no_trials)=100; %There are no trials before the first trial
                    else
                        handles.drgb.PAC.delta_t_trial(no_trials)=handles.drg.session(sessionNo).trial_start(trNo)-handles.drg.session(sessionNo).trial_start(trNo-1);
                    end
                    perCorr_per_histo(no_trials)=50;
                    
                    if handles.peakLFPNo==18
                        %This is the sniff
                        [meanVectorLength(no_trials), meanVectorAngle(no_trials), peakAngle(no_trials), mod_indx(no_trials), phase, phase_histo, theta_wave]=drgGetThetaAmpPhaseSniff(LFPlow,LFPhigh,Fs,lowF1,lowF2,highF1,highF2,pad_time,n_phase_bins,handles.which_method);
                    else
                        [meanVectorLength(no_trials), meanVectorAngle(no_trials), peakAngle(no_trials), mod_indx(no_trials), phase,...
                            phase_histo, theta_wave, meanPeakAngle, out_times, out_phase, out_time_PAChisto, decLFPgenv, decanglethetaLFP, out_times_env]...
                            =drgGetThetaAmpPhase(LFPlow,LFPhigh,Fs,lowF1,lowF2,highF1,highF2,pad_time,n_phase_bins,handles.which_method);
                    end
                    
                    conc_evTyNoPerTr(no_trials)=conc_evTyNo;
                    
                    out_times=out_times+handles.time_start+handles.time_pad;
                    
                    %Save the output
                    handles.drgb.PAC.no_trials=no_trials;
                    handles.drgb.PAC.meanVectorLength(no_trials)=meanVectorLength(no_trials);
                    handles.drgb.PAC.meanVectorAngle(no_trials)=meanVectorAngle(no_trials);
                    handles.drgb.PAC.peakAngle(no_trials)=peakAngle(no_trials);
                    handles.drgb.PAC.mod_indx(no_trials)=mod_indx(no_trials);
                    handles.drgb.PAC.all_phase_histo(no_trials,1:n_phase_bins+1)=phase_histo;
                    handles.drgb.PAC.all_theta_wave(no_trials,1:n_phase_bins+1)=theta_wave;
                    handles.drgb.PAC.perCorr(no_trials)=perCorr(drgFindEvNo(handles,trNo,sessionNo,odorOn));
                    handles.drgb.PAC.meanPeakAngle(no_trials)=meanPeakAngle;
                    handles.drgb.PAC.PACtimecourse(no_trials).out_times=out_times;
                    handles.drgb.PAC.PACtimecourse(no_trials).out_phase=out_phase;
                    handles.drgb.PAC.PACtimecourse(no_trials).out_time_PAChisto=out_time_PAChisto;
                    handles.drgb.PAC.PACtimecourse(no_trials).decLFPgenv=decLFPgenv;
                    handles.drgb.PAC.PACtimecourse(no_trials).decanglethetaLFP=decanglethetaLFP;
                    handles.drgb.PAC.out_times_env=out_times_env;
                    handles.drgb.PAC.phase=phase;
                    
                    all_out_phase(no_trials,1:length(out_times))=out_phase;
                    if no_trials==1
                        all_out_time_PAChisto=out_time_PAChisto;
                    else
                        all_out_time_PAChisto=all_out_time_PAChisto+out_time_PAChisto;
                    end
                    all_phase_histo(no_trials,1:n_phase_bins+1)=phase_histo;
                    all_theta_wave(no_trials,1:n_phase_bins+1)=theta_wave;
                    no_encoding_trials=no_encoding_trials+1;
                    enc_phase_histo(no_encoding_trials,1:n_phase_bins+1)=phase_histo;
                    MI_enc=[MI_enc mod_indx(no_trials)];
                    
                    if handles.save_drgb==1
                        for evTypeNo=1:length(handles.drgbchoices.evTypeNos)
                            if sum(handles.drg.session(1).events(handles.drgbchoices.evTypeNos(evTypeNo)).times==handles.drg.session(1).events(handles.drgbchoices.referenceEvent).times(evNo))>0
                                handles.drgb.PAC.which_event(evTypeNo,no_trials)=1;
                            else
                                handles.drgb.PAC.which_event(evTypeNo,no_trials)=0;
                            end
                        end
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
    
end

handles.evTypeNo=currentevTypeNo;

if ~isempty(spm)
    %Now figure out the gamma power at the peak
    if handles.displayData==0
        if handles.use_peakAngle==0
            %With no display (i.e., for batch) do this calculation only for S+
            %Note that I use phase, mean(enc_phase_histo), as opposed to
            %meanPeakAngle because meanPeakAngle may switch spuriously from
            %trial to trial when there are two peaks
            
            [max_hist max_ii]=max(mean(enc_phase_histo(logical(spm),:),1));
            meanPeakAngle=(phase(max_ii)*(pi/180))-pi;
            handles.drgb.PAC.peakAngleForPower=phase(max_ii);
            
            [min_hist min_ii]=min(mean(enc_phase_histo(logical(spm),:),1));
            meanTroughAngle=(phase(min_ii)*(pi/180))-pi;
            handles.drgb.PAC.troughAngleForPower=phase(min_ii);
        else
            meanPeakAngle=(handles.peakAngle_for_power*(pi/180))-pi;
            handles.drgb.PAC.peakAngleForPower=handles.peakAngle_for_power;
            meanTroughAngle=(handles.troughAngle_for_power*(pi/180))-pi;
            handles.drgb.PAC.troughAngleForPower=handles.troughAngle_for_power;
        end
    else
        if handles.use_peakAngle==0
            %Note that I use phase, mean(enc_phase_histo), as opposed to
            %meanPeakAngle because meanPeakAngle may switch spuriously from
            %trial to trial when there are two peaks
            [max_hist max_ii]=max(mean(enc_phase_histo));
            meanPeakAngle=(phase(max_ii)*(pi/180))-pi;
            handles.drgb.PAC.peakAngleForPower=phase(max_ii);
            if handles.displayData==1
                handles.peakAngle_for_power=phase(max_ii);
                set(handles.set_peakAngle,'String',num2str(phase(max_ii)));
            end
            [min_hist min_ii]=min(mean(enc_phase_histo));
            meanTroughAngle=(phase(min_ii)*(pi/180))-pi;
            handles.drgb.PAC.troughAngleForPower=phase(min_ii);
            if handles.displayData==1
                handles.troughAngle_for_power=phase(min_ii);
                set(handles.set_troughAngle,'String',num2str(phase(min_ii)));
            end
        else
            meanPeakAngle=(handles.peakAngle_for_power*(pi/180))-pi;
            handles.drgb.PAC.peakAngleForPower=handles.peakAngle_for_power;
            meanTroughAngle=(handles.troughAngle_for_power*(pi/180))-pi;
            handles.drgb.PAC.troughAngleForPower=handles.troughAngle_for_power;
        end
    end
    
    dt=handles.window-handles.noverlap;
    t=[handles.time_start+handles.time_pad:dt:handles.time_end-handles.time_pad];
    handles.drgb.PAC.t_power=t;
    peakPower=zeros(no_trials,length(t));
    %Find the value of gamma power at each point
    out_times_env=out_times_env+handles.time_start+handles.time_pad;
    
    %Find peak power
    for trNum=1:no_trials
        peakAmp=[];
        peakAmp_times=[];
        ii=1;
        jj=0;
        at_end=0;
        this_angleTetaLFP=handles.drgb.PAC.PACtimecourse(trNum).decanglethetaLFP;
        this_LFPenv=handles.drgb.PAC.PACtimecourse(trNum).decLFPgenv;
        while at_end==0
            ii_next=find(this_angleTetaLFP(ii:end)>=meanPeakAngle,1,'first');
            if (~isempty(ii_next))&(ii+ii_next-1<=length(out_times_env))
                jj=jj+1;
                peakAmp(jj)=this_LFPenv(ii+ii_next-1);
                peakAmp_times(jj)=out_times_env(ii+ii_next-1);
                ii=ii+ii_next;
                ii_next=find(this_angleTetaLFP(ii:end)<meanPeakAngle,1,'first');
                if ~isempty(ii_next)
                    ii=ii+ii_next;
                else
                    at_end=1;
                end
            else
                at_end=1;
            end
        end
        
        for ii_t=1:length(t)
            if t(ii_t)<=peakAmp_times(1)
                peakPower(trNum,ii_t)=10*log10((peakAmp(1)^2)/2);
                %peakAmp_t(ii_t)=peakAmp(1);
            else
                if t(ii_t)>=peakAmp_times(end)
                    peakPower(trNum,ii_t)=10*log10((peakAmp(end)^2)/2);
                    %peakAmp_t(ii_t)=peakAmp(end);
                else
                    ii_pt=find(peakAmp_times>=t(ii_t),1,'first');
                    this_amp=peakAmp(ii_pt-1)+(t(ii_t)-peakAmp_times(ii_pt-1))*((peakAmp(ii_pt)-peakAmp(ii_pt-1))/(peakAmp_times(ii_pt)-peakAmp_times(ii_pt-1)));
                    %peakAmp_t(ii_t)=this_amp;
                    peakPower(trNum,ii_t)=10*log10((this_amp^2)/2);
                end
            end
        end
        
        if handles.subtractRef==1
            peakPower(trNum,:)=peakPower(trNum,:)-mean(peakPower(trNum,(t>=handles.startRef+handles.time_pad)&(t<=handles.endRef-handles.time_pad)));
        end
        
        handles.drgb.PAC.PACtimecourse(trNum).peakPower=peakPower(trNum,:);
        handles.drgb.PAC.meanPeakPower(trNum)=mean(peakPower(trNum,:));
    end
    
    %Find trough power
    troughPower=zeros(no_trials,length(t));
    for trNum=1:no_trials
        troughAmp=[];
        troughAmp_times=[];
        ii=1;
        jj=0;
        at_end=0;
        this_angleTetaLFP=handles.drgb.PAC.PACtimecourse(trNum).decanglethetaLFP;
        this_LFPenv=handles.drgb.PAC.PACtimecourse(trNum).decLFPgenv;
        while at_end==0
            ii_next=find(this_angleTetaLFP(ii:end)>=meanTroughAngle,1,'first');
            if (~isempty(ii_next))&(ii+ii_next-1<=length(out_times_env))
                jj=jj+1;
                troughAmp(jj)=this_LFPenv(ii+ii_next-1);
                troughAmp_times(jj)=out_times_env(ii+ii_next-1);
                ii=ii+ii_next;
                ii_next=find(this_angleTetaLFP(ii:end)<meanTroughAngle,1,'first');
                if ~isempty(ii_next)
                    ii=ii+ii_next;
                else
                    at_end=1;
                end
            else
                at_end=1;
            end
        end
        
        for ii_t=1:length(t)
            if t(ii_t)<=troughAmp_times(1)
                troughPower(trNum,ii_t)=10*log10((troughAmp(1)^2)/2);
                %troughAmp_t(ii_t)=troughAmp(1);
            else
                if t(ii_t)>=troughAmp_times(end)
                    troughPower(trNum,ii_t)=10*log10((troughAmp(end)^2)/2);
                    %troughAmp_t(ii_t)=troughAmp(end);
                else
                    ii_pt=find(troughAmp_times>=t(ii_t),1,'first');
                    this_amp=troughAmp(ii_pt-1)+(t(ii_t)-troughAmp_times(ii_pt-1))*((troughAmp(ii_pt)-troughAmp(ii_pt-1))/(troughAmp_times(ii_pt)-troughAmp_times(ii_pt-1)));
                    %troughAmp_t(ii_t)=this_amp;
                    troughPower(trNum,ii_t)=10*log10((this_amp^2)/2);
                end
            end
        end
        
        if handles.subtractRef==1
            troughPower(trNum,:)=troughPower(trNum,:)-mean(troughPower(trNum,(t>=handles.startRef+handles.time_pad)&(t<=handles.endRef-handles.time_pad)));
        end
        
        handles.drgb.PAC.PACtimecourse(trNum).troughPower=troughPower(trNum,:);
        handles.drgb.PAC.meanTroughPower(trNum)=mean(troughPower(trNum,:));
    end
    
    
    if handles.displayData==1
        
        all_out_time_PAChisto=all_out_time_PAChisto/no_trials;
        mean_MI_enc=mean(MI_enc);
        
          conc_labels{16}='10';
        conc_labels{17}='3.2';
        conc_labels{18}='1';
        conc_labels{19}='0.32';
        conc_labels{20}='0.1';
        conc_labels{21}='0.032';
        
        %Plot the PAC per trial for S+
        try
            close 1
        catch
        end
        
        hFig1 = figure(1);
        set(hFig1, 'units','normalized','position',[.01 .1 .23 .8])
        
        min_prob=prctile(all_phase_histo(:),5);
        max_prob=prctile(all_phase_histo(:),95);
        
        splus_trials=(conc_evTyNoPerTr==16)|(conc_evTyNoPerTr==17)|(conc_evTyNoPerTr==18);
        nosp=sum(splus_trials);
        
        trials=1:nosp;
        %pcolor(repmat(phase,length(trials),1),repmat(trials',1,length(phase)),all_phase_histo)
        drg_pcolor(repmat(phase,length(trials),1),repmat(trials',1,length(phase)),all_phase_histo(splus_trials,:))
        colormap jet
        shading flat
        % min_prob=0.0113;
        % max_prob=0.0314;
        caxis([min_prob    max_prob])
        xlabel('Phase for low freq oscillation (deg)')
        ylabel('Trial No');
        title(['Phase-amplitude coupling per trial S+'])
        
        for conc_evNo=16:18
            firstTr=find(conc_evTyNoPerTr==conc_evNo,1,'first')-find(conc_evTyNoPerTr==18,1,'first')+1;
            lastTr=find(conc_evTyNoPerTr==conc_evNo,1,'last')-find(conc_evTyNoPerTr==18,1,'first')+1;
            fprintf(1, ['\nFor ' conc_labels{conc_evNo} '%% from trial No %d to %d\n'], firstTr,lastTr);
        end
        
        fprintf(1, ['\n\n']);
        
        %Plot the PAC per trial for S-
        try
            close 2
        catch
        end
        
        hFig2 = figure(2);
        set(hFig2, 'units','normalized','position',[.01 .1 .23 .8])
        
        sminus_trials=(conc_evTyNoPerTr==19)|(conc_evTyNoPerTr==20)|(conc_evTyNoPerTr==21);
        nosm=sum(sminus_trials);
        
        trials=1:nosm;
        %pcolor(repmat(phase,length(trials),1),repmat(trials',1,length(phase)),all_phase_histo)
        drg_pcolor(repmat(phase,length(trials),1),repmat(trials',1,length(phase)),all_phase_histo(sminus_trials,:))
        colormap jet
        shading flat
        % min_prob=0.0113;
        % max_prob=0.0314;
        caxis([min_prob    max_prob])
        xlabel('Phase for low freq oscillation (deg)')
        ylabel('Trial No');
        title(['Phase-amplitude coupling per trial S-'])
        
           for conc_evNo=19:21
            firstTr=find(conc_evTyNoPerTr==conc_evNo,1,'first');
            lastTr=find(conc_evTyNoPerTr==conc_evNo,1,'last');
            fprintf(1, ['\nFor ' conc_labels{conc_evNo} '%% from trial No %d to %d\n'], firstTr,lastTr);
        end
        
        fprintf(1, ['\n\n']);
        
        try
            close 3
        catch
        end
        
        hFig3 = figure(3);
        set(hFig3, 'units','normalized','position',[.49 .1 .05 .3])
        
        
        prain=[min_prob:(max_prob-min_prob)/99:max_prob];
        pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
        colormap jet
        shading interp
        ax=gca;
        set(ax,'XTickLabel','')
        
        %Plot a bar graph for MI per conc event type
        try
            close 4
        catch
        end
        
        hFig4 = figure(4);
        set(hFig4, 'units','normalized','position',[.15 .6 .3 .3])
        
        hold on
        
        
        position=0;
        
        for conc_evTyNo=last_conc_evNo:-1:first_conc_evNo
            position=position+1;
            mean_MI=mean(mod_indx((conc_evTyNoPerTr==conc_evTyNo)));
            bar(position,mean_MI,'FaceColor',[0.7 0.7 1])
            CI_MI = bootci(1000, {@mean, mod_indx((conc_evTyNoPerTr==conc_evTyNo))})';
            plot([position position],CI_MI,'-k','LineWidth',3)
            plot(position*ones(1,sum(conc_evTyNoPerTr==conc_evTyNo)),mod_indx((conc_evTyNoPerTr==conc_evTyNo)),'ok')
        end
        
        ylabel('MI')
        xticks([1:6])
        xticklabels({'0.032','0.1','0.32','1','3.2','10'})
        xlabel('Percent isoamyl acetate')
        
        glm_ii=0;
        glm_mi=[];
        ii=0;
        input_data=[];
        
    
        
        for conc_evTyNo=last_conc_evNo:-1:first_conc_evNo
            %Save data for glm and ranksum
            glm_mi.data(glm_ii+1:glm_ii+sum((conc_evTyNoPerTr==conc_evTyNo)))=mod_indx((conc_evTyNoPerTr==conc_evTyNo));
            if conc_evTyNo<=18
                glm_mi.spm(glm_ii+1:glm_ii+sum((conc_evTyNoPerTr==conc_evTyNo)))=zeros(1,sum((conc_evTyNoPerTr==conc_evTyNo)));
            else
                glm_mi.spm(glm_ii+1:glm_ii+sum((conc_evTyNoPerTr==conc_evTyNo)))=ones(1,sum((conc_evTyNoPerTr==conc_evTyNo)));
            end
            glm_mi.conc(glm_ii+1:glm_ii+sum((conc_evTyNoPerTr==conc_evTyNo)))=conc_evTyNo*ones(1,sum((conc_evTyNoPerTr==conc_evTyNo)));
            glm_ii=glm_ii+sum((conc_evTyNoPerTr==conc_evTyNo));
            
            ii=ii+1;
            input_data(ii).data=mod_indx((conc_evTyNoPerTr==conc_evTyNo));
            input_data(ii).description=[conc_labels{conc_evTyNo}];
        end
        
        %Perform the glm
        fprintf(1, ['glm for MI\n'])
        
        
        %There is only one group here (e.g. for Justin's paper we only include
        %forward)
        fprintf(1, ['\n\nglm for MI\n'])
        tbl = table(glm_mi.data',glm_mi.spm',glm_mi.conc',...
            'VariableNames',{'MI','spm','conc'});
        mdl = fitglm(tbl,'MI~spm+conc'...
            ,'CategoricalVars',[2,3])
        
        
        %Do the ranksum/t-test
        fprintf(1, ['\n\nRanksum or t-test p values for  MI\n'])
        [output_data] = drgMutiRanksumorTtest(input_data);
        
        
        %Plot the peak angle circular histos
        
        try
            close 5
        catch
        end
        
        hFig5 = figure(5);
        set(hFig5, 'units','normalized','position',[.39 .05 .45 .45])
        
        for conc_evTyNo=last_conc_evNo:-1:first_conc_evNo
            subplot(2,3,conc_evTyNo-15)
            
            polarhistogram(pi*peakAngle((conc_evTyNoPerTr==conc_evTyNo))/180,12)
            title(['Peak angle ' conc_labels{conc_evTyNo}])
            
        end
        
        
        pffft=1;
    end
end

if handles.displayData==1
    fprintf(1, ['\nPCA processed for %d out of %d trials \n\n'], no_trials,trials_attempted);
end



pfffft=1;






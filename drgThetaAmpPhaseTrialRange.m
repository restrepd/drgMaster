function handles=drgThetaAmpPhaseTrialRange(handles)

%Generates a trial per trial phase histogram
odorOn=2;
splus=5;
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


for trNo=firstTr:lastTr
    
    if handles.displayData==1
        trial_no=trNo
    end
    
    evNo = drgFindEvNo(handles,trNo,sessionNo);
    
    if evNo~=-1
        excludeTrial=drgExcludeTrialLFP(handles.drg,handles.peakLFPNo,handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);
        
        if excludeTrial==0
            
            [LFPlow, trialNo, can_read1] = drgGetTrialLFPData(handles, handles.peakLFPNo, evNo, handles.evTypeNo, handles.time_start, handles.time_end);
            [LFPhigh, trialNo, can_read2] = drgGetTrialLFPData(handles, handles.burstLFPNo, evNo, handles.evTypeNo, handles.time_start, handles.time_end);
            
            if (can_read1==1)&(can_read2==1)
                no_trials=no_trials+1;
               %If this is being used for batch processing find out whether
                %this is an S+
                eventNo = drgFindEvNo(handles,trNo,sessionNo,splus);
                if eventNo~=-1
                    spm(no_trials)=1;
                else
                    spm(no_trials)=0;
                end
                
                time(no_trials)=handles.drg.session(sessionNo).trial_start(trialNo);
                handles.drgb.PAC.this_trialNo(no_trials)=trNo;
                perCorr_per_histo(no_trials)=50;
                if handles.peakLFPNo==18
                    %This is the sniff
                    [meanVectorLength(no_trials), meanVectorAngle(no_trials), peakAngle(no_trials), mod_indx(no_trials), phase, phase_histo, theta_wave]=drgGetThetaAmpPhaseSniff(LFPlow,LFPhigh,Fs,lowF1,lowF2,highF1,highF2,pad_time,n_phase_bins,handles.which_method);
                else
                    [meanVectorLength(no_trials), meanVectorAngle(no_trials), peakAngle(no_trials), mod_indx(no_trials), phase,...
                        phase_histo, theta_wave, meanPeakAngle, out_times, out_phase, out_time_PAChisto, decLFPgenv, decanglethetaLFP, out_times_env]...
                        =drgGetThetaAmpPhase(LFPlow,LFPhigh,Fs,lowF1,lowF2,highF1,highF2,pad_time,n_phase_bins,handles.which_method);
                end
                
                out_times=out_times+handles.time_start+handles.time_pad;
                
                %Save the output
                handles.drgb.PAC.no_trials=no_trials;
                handles.drgb.PAC.meanVectorLength(no_trials)=meanVectorLength(no_trials);
                handles.drgb.PAC.meanVectorAngle(no_trials)=meanVectorAngle(no_trials);
                handles.drgb.PAC.peakAngle(no_trials)=peakAngle(no_trials);
                handles.drgb.PAC.mod_indx(no_trials)=mod_indx(no_trials);
                handles.drgb.PAC.all_phase_histo(no_trials,1:n_phase_bins+1)=phase_histo;
                handles.drgb.PAC.all_theta_wave(no_trials,1:n_phase_bins+1)=theta_wave;
                handles.drgb.PAC.perCorr(no_trials)=perCorr(drgFindEvNo(handles,trialNo,sessionNo,odorOn));
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
            end
        end
    end
    %end
    %end %if eventstamps...
end %for evNo

%Now figure out the gamma power at the peak
if handles.displayData==0
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
    if handles.use_peakAngle==0
        %Note that I use phase, mean(enc_phase_histo), as opposed to
        %meanPeakAngle because meanPeakAngle may switch spuriously from
        %trial to trial when there are two peaks
        [max_hist max_ii]=max(mean(enc_phase_histo));
        meanPeakAngle=(phase(max_ii)*(pi/180))-pi;
        if handles.displayData==1
            handles.peakAngle_for_power=phase(max_ii);
            set(handles.set_peakAngle,'String',num2str(phase(max_ii)));
        end
        [min_hist min_ii]=min(mean(enc_phase_histo));
        meanTroughAngle=(phase(min_ii)*(pi/180))-pi;
        if handles.displayData==1
            handles.troughAngle_for_power=phase(min_ii);
            set(handles.set_troughAngle,'String',num2str(phase(min_ii)));
        end
    else
        meanPeakAngle=(handles.peakAngle_for_power*(pi/180))-pi;
        meanTroughAngle=(handles.troughAngle_for_power*(pi/180))-pi;
    end
end

dt=handles.window-handles.noverlap;
t=[handles.time_start+handles.time_pad:dt:handles.time_end-handles.time_pad];
handles.drgb.PAC.t_power=t;
peakPower=zeros(no_trials,length(t));
%Find the value of gamma power at each point
out_times_env=out_times_env+handles.time_start+handles.time_pad;

%Find peak power
for trNo=1:no_trials
    peakAmp=[];
    peakAmp_times=[];
    ii=1;
    jj=0;
    at_end=0;
    this_angleTetaLFP=handles.drgb.PAC.PACtimecourse(trNo).decanglethetaLFP;
    this_LFPenv=handles.drgb.PAC.PACtimecourse(trNo).decLFPgenv;
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
            peakPower(trNo,ii_t)=10*log10((peakAmp(1)^2)/2);
            %peakAmp_t(ii_t)=peakAmp(1);
        else
            if t(ii_t)>=peakAmp_times(end)
                peakPower(trNo,ii_t)=10*log10((peakAmp(end)^2)/2);
                %peakAmp_t(ii_t)=peakAmp(end);
            else
                ii_pt=find(peakAmp_times>=t(ii_t),1,'first');
                this_amp=peakAmp(ii_pt-1)+(t(ii_t)-peakAmp_times(ii_pt-1))*((peakAmp(ii_pt)-peakAmp(ii_pt-1))/(peakAmp_times(ii_pt)-peakAmp_times(ii_pt-1)));
                %peakAmp_t(ii_t)=this_amp;
                peakPower(trNo,ii_t)=10*log10((this_amp^2)/2);
            end
        end
    end
    
    if handles.subtractRef==1
        peakPower(trNo,:)=peakPower(trNo,:)-mean(peakPower(trNo,(t>=handles.startRef+handles.time_pad)&(t<=handles.endRef-handles.time_pad)));
    end
    
    handles.drgb.PAC.PACtimecourse(trNo).peakPower=peakPower(trNo,:);
    handles.drgb.PAC.meanPeakPower(trNo)=mean(peakPower(trNo,:));
end

%Find trough power
troughPower=zeros(no_trials,length(t));
for trNo=1:no_trials
    troughAmp=[];
    troughAmp_times=[];
    ii=1;
    jj=0;
    at_end=0;
    this_angleTetaLFP=handles.drgb.PAC.PACtimecourse(trNo).decanglethetaLFP;
    this_LFPenv=handles.drgb.PAC.PACtimecourse(trNo).decLFPgenv;
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
            troughPower(trNo,ii_t)=10*log10((troughAmp(1)^2)/2);
            %troughAmp_t(ii_t)=troughAmp(1);
        else
            if t(ii_t)>=troughAmp_times(end)
                troughPower(trNo,ii_t)=10*log10((troughAmp(end)^2)/2);
                %troughAmp_t(ii_t)=troughAmp(end);
            else
                ii_pt=find(troughAmp_times>=t(ii_t),1,'first');
                this_amp=troughAmp(ii_pt-1)+(t(ii_t)-troughAmp_times(ii_pt-1))*((troughAmp(ii_pt)-troughAmp(ii_pt-1))/(troughAmp_times(ii_pt)-troughAmp_times(ii_pt-1)));
                %troughAmp_t(ii_t)=this_amp;
                troughPower(trNo,ii_t)=10*log10((this_amp^2)/2);
            end
        end
    end
    
    if handles.subtractRef==1
        troughPower(trNo,:)=troughPower(trNo,:)-mean(troughPower(trNo,(t>=handles.startRef+handles.time_pad)&(t<=handles.endRef-handles.time_pad)));
    end
    
    handles.drgb.PAC.PACtimecourse(trNo).troughPower=troughPower(trNo,:);
    handles.drgb.PAC.meanTroughPower(trNo)=mean(troughPower(trNo,:));
end


if handles.displayData==1
    
    all_out_time_PAChisto=all_out_time_PAChisto/no_trials; 
    mean_MI_enc=mean(MI_enc);
    
    %Plot the timecourse for the phase
    try
        close 1
    catch
    end
    
    hFig1 = figure(1);
    set(hFig1, 'units','normalized','position',[.15 .6 .8 .23])
    
    %     plot(out_times,circ_rad2ang(circ_mean(circ_ang2rad(all_out_phase)-180))+180,'ob')
    plot(out_times,(180/pi)*circ_axial(circ_mean(all_out_phase*pi/180)'),'ob')
    ylim([0 360])
    
    
    mean_peakPower=mean(peakPower,1)';
    CI_peakPower = bootci(1000, {@mean, peakPower})';
    CI_peakPower(:,1)=mean_peakPower-CI_peakPower(:,1);
    CI_peakPower(:,2)=CI_peakPower(:,2)-mean_peakPower;
    

    mean_troughPower=mean(troughPower,1)';
    CI_troughPower = bootci(1000, {@mean, troughPower})';
    CI_troughPower(:,1)=mean_troughPower-CI_troughPower(:,1);
    CI_troughPower(:,2)=CI_troughPower(:,2)-mean_troughPower;
    
    try
        close 7
    catch
    end
    
    hFig7 = figure(7);
    set(hFig7, 'units','normalized','position',[.25 .2 .7 .23])
    hold on
    
    [hltrough, hptrough]=boundedline(t,mean_troughPower, CI_troughPower, 'b');
    [hlpeak, hppeak]=boundedline(t,mean_peakPower, CI_peakPower, 'r');
    xlabel('Time(sec)')
    ylabel('Peak power (dB)')
    legend([hltrough hlpeak],{'Trough','Peak'})
    
    if handles.autoscale==0
        ylim([handles.minLogP handles.maxLogP]);
    end
    
    try
        close 6
    catch
    end
    
    hFig6 = figure(6);
    set(hFig6, 'units','normalized','position',[.25 .2 .7 .23])
    
    min_prob=prctile(all_out_time_PAChisto(:),5);
    max_prob=prctile(all_out_time_PAChisto(:),95);
    
    
    %pcolor(repmat(phase,length(trials),1),repmat(trials',1,length(phase)),all_phase_histo)
    drg_pcolor(repmat(out_times,length(phase(1:end-1)),1),repmat(phase(1:end-1)',1,length(out_times)),all_out_time_PAChisto')
    colormap jet
    shading flat
    % min_prob=0.0113;
    % max_prob=0.0314;
    caxis([min_prob    max_prob])
    xlabel('Time (sec)')
    ylabel('Phase for low freq oscillation (deg)');
    title(['Phase-amplitude coupling timecourse ' handles.drg.session(1).draq_d.eventlabels{handles.evTypeNo}])
    
    %     (180/pi)*circ_mean(all_out_phase'*pi/180)'
    try
        close 2
    catch
    end
    
    hFig2 = figure(2);
    set(hFig2, 'units','normalized','position',[.25 .1 .23 .8])
    
    %Plot the modulation index as a function of trial number
    subplot(4,1,1)
    trials=1:no_trials;
    plot(trials,mod_indx,'o-')
    hold on
    if length(mod_indx>20)
        no_conv_points=6;
        conv_win=ones(1,no_conv_points);
        mod_indx_extend=[mean(mod_indx(1:no_conv_points/2))*ones(1,no_conv_points) mod_indx mean(mod_indx(end-(no_conv_points/2)+1:end))*ones(1,no_conv_points)];
        conv_mod_indx_extend=conv(mod_indx_extend,conv_win)/no_conv_points;
        conv_mod_indx=conv_mod_indx_extend(no_conv_points+1:no_conv_points+length(mod_indx));
        plot(trials,conv_mod_indx,'-b','LineWidth',3)
    end
    xlabel('Trial No')
    ylabel('Modulation index')
    title('Modulation index vs trial number')
    
    % figure(11)
    % plot(trials,mod_indx,'o-k','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',7,'LineWidth',2)
    % xlabel('Trial No')
    % ylabel('Modulation index')
    % title('Modulation index vs trial number')
    % ylim([0 0.07])
    % set(gca,'FontSize',20,'LineWidth',1.25)
    
    
    
    %Now plot the mean theta waveform
    subplot(4,1,2)
    shadedErrorBar(phase,mean(all_theta_wave,1),std(all_theta_wave,0,1),'-b')
    xlim([0 360])
    title('Mean low frequency waveform')
    xlabel('Degrees')
    
    %Now plot the encoding theta phase histogram for gamma
    subplot(4,1,3)
    shadedErrorBar(phase,mean(enc_phase_histo,1),std(enc_phase_histo,0,1),'-b')
    xlim([0 360])
    title('Phase-amplitude coupling')
    ylabel('Probability')
    xlabel('Phase for low freq oscillation (deg)')
    yl_enc=ylim;
    text(20,yl_enc(1)+0.9*(yl_enc(2)-yl_enc(1)),['MI= ' num2str(mean_MI_enc)])
    
    
    
    try
        close 3
    catch
    end
    
    hFig3 = figure(3);
    set(hFig3, 'units','normalized','position',[.01 .1 .23 .8])
    
    min_prob=prctile(all_phase_histo(:),5);
    max_prob=prctile(all_phase_histo(:),95);
    
    trials=1:no_trials;
    %pcolor(repmat(phase,length(trials),1),repmat(trials',1,length(phase)),all_phase_histo)
    drg_pcolor(repmat(phase,length(trials),1),repmat(trials',1,length(phase)),all_phase_histo)
    colormap jet
    shading flat
    % min_prob=0.0113;
    % max_prob=0.0314;
    caxis([min_prob    max_prob])
    xlabel('Phase for low freq oscillation (deg)')
    ylabel('Trial No');
    title(['Phase-amplitude coupling per trial' handles.drg.session(1).draq_d.eventlabels{handles.evTypeNo}])
    
    try
        close 4
    catch
    end
    
    hFig4 = figure(4);
    set(hFig4, 'units','normalized','position',[.49 .1 .05 .3])
    
    
    prain=[min_prob:(max_prob-min_prob)/99:max_prob];
    pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
    colormap jet
    shading interp
    ax=gca;
    set(ax,'XTickLabel','')
    
    try
        close 5
    catch
    end
    
    hFig5 = figure(5);
    set(hFig5, 'units','normalized','position',[.59 .05 .35 .35])
    
    if no_encoding_trials>0
        rose(pi*meanVectorAngle/180,12)
        title('Mean phase')
    end
    
    
    
    pffft=1;
end

if handles.displayData==1
    fprintf(1, ['\nPCA processed for %d out of %d trials \n\n'], no_trials,lastTr);
end

pfffft=1;






function handles=drgThetaAmpPhaseLick(handles)

%Generates a trial per trial phase histogram
sessionNo=handles.sessionNo;
Fs=handles.drg.session(sessionNo).draq_p.ActualRate;
lowF1=handles.peakLowF;
lowF2=handles.peakHighF;
highF1=handles.burstLowF;
highF2=handles.burstHighF;
pad_time=handles.time_pad;
n_phase_bins=handles.n_phase_bins;

%Enter trials
firstTr=handles.trialNo;
lastTr=handles.lastTrialNo;
odorOn=2;

[perCorr, encoding_trials, retrieval_trials, encoding_this_evTypeNo,retrieval_this_evTypeNo]=drgFindEncRetr(handles);

try
    close 1
catch
end

hFig1 = figure(1);
set(hFig1, 'units','normalized','position',[.55 .27 .35 .15])


% %Plot the percent correct
% trials=1:length(perCorr);
% plot(trials,perCorr,'o','MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7])
% hold on
% plot(trials(encoding_trials),perCorr(encoding_trials),'ob')
% plot(trials(retrieval_trials),perCorr(retrieval_trials),'or')
% ylim([40 110]);
% xlabel('Trial No')
% ylabel('Percent correct')
% title('% correct vs trial number, blue=encoding, red=retrieval')

%Now get the phase of gamma bursts witin theta
no_encoding_trials=0;
no_retrieval_trials=0;
no_trials=0;
enc_phase_histo=[];
retr_phase_histo=[];
all_phase_histo=[];
all_theta_wave=[];
MI_enc=[];
MI_retr=[];
which_event=[];
graphHit=[];
graphMiss=[];
graphCR=[];
graphFA=[];

handles=drgPercentLickPerTrial(handles);
handles=drgBusquetAnalysis(handles);

for trNo=firstTr:lastTr
    
    if handles.save_drgb==0
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
                time(no_trials)=handles.drg.session(sessionNo).trial_start(trialNo);
                which_trial(no_trials)=1;
                perCorr_per_histo(no_trials)= perCorr(drgFindEvNo(handles,trialNo,sessionNo,odorOn));
                [meanVectorLength(no_trials), meanVectorAngle(no_trials), peakAngle(no_trials), mod_indx(no_trials), phase, phase_histo, theta_wave]=drgGetThetaAmpPhase(LFPlow,LFPhigh,Fs,lowF1,lowF2,highF1,highF2,pad_time,n_phase_bins,handles.which_method);
                all_phase_histo(no_trials,1:n_phase_bins+1)=phase_histo;
                all_theta_wave(no_trials,1:n_phase_bins+1)=theta_wave;
                no_encoding_trials=no_encoding_trials+1;
                enc_phase_histo(no_encoding_trials,1:n_phase_bins+1)=phase_histo;
                MI_enc=[MI_enc mod_indx(no_trials)];
                percent_lick(no_trials)=handles.drg.session(sessionNo).percent_lick(trNo);
                graphHit(no_trials)=handles.gonogo.graphHit(trNo);
                graphMiss(no_trials)=handles.gonogo.graphMiss(trNo);
                graphCR(no_trials)=handles.gonogo.graphCR(trNo);
                graphFA(no_trials)=handles.gonogo.graphFA(trNo);
                
                if handles.save_events==1
                    for evTypeNo=1:length(handles.drgb.choices.evTypeNos)
                        if sum(handles.drg.session(1).events(handles.drgb.choices.evTypeNos(evTypeNo)).times==handles.drg.session(1).events(handles.drgb.choices.referenceEvent).times(evNo))>0
                            which_event(evTypeNo,no_trials)=1;
                        else
                            which_event(evTypeNo,no_trials)=0;
                        end
                    end
                end
            end
        end
    end
    %end
    %end %if eventstamps...
end %for evNo

mean_MI_enc=mean(MI_enc);


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
xlabel('Trial No')
ylabel('Modulation index')
title('Modulation index vs trial number')

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



trials=1:no_trials;
%pcolor(repmat(phase,length(trials),1),repmat(trials',1,length(phase)),all_phase_histo)
drg_pcolor(repmat(phase,length(trials),1),repmat(trials',1,length(phase)),all_phase_histo)
colormap jet
shading flat
z_per=prctile(all_phase_histo(:),99);
z_min=prctile(all_phase_histo(:),1);
caxis([z_min z_per]);
xlabel('Phase for low freq oscillation (deg)')
ylabel('Trial No');
title(['Phase-amplitude coupling per trial' handles.drg.session(1).draq_d.eventlabels{handles.evTypeNo}])

try
    close 4
catch
end

hFig4 = figure(4);
set(hFig4, 'units','normalized','position',[.49 .1 .05 .3])

min_prob=min(min(all_phase_histo));
max_prob=max(max(all_phase_histo));
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
set(hFig5, 'units','normalized','position',[.49 .05 .35 .35])

if no_encoding_trials>0
    rose(pi*meanVectorAngle/180,12)
    title('Mean phase')
end

%Save the data if requested by drgb
if handles.save_drgb==1
    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no).no_trials=no_trials;
    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no).which_trial=which_trial;
    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no).meanVectorLength=meanVectorLength;
    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no).meanVectorAngle=meanVectorAngle;
    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no).peakAngle=peakAngle';
    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no).mod_indx=mod_indx;
    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no).all_phase_histo=all_phase_histo;
    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no).perCorr_per_histo=perCorr_per_histo;
    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no).which_event=which_event;
    handles.drgb.phase=phase;
end



try
    close 6
catch
end

hFig6 = figure(6);
set(hFig6, 'units','normalized','position',[.49 .05 .35 .35])



legendInfo=[];
li_ii=0;
p1=polarplot(pi*meanVectorAngle(graphHit==1)/180,percent_lick(graphHit==1),' vr', 'MarkerFaceColor','r', 'MarkerSize',7)
if ~isempty(p1)
    li_ii=li_ii+1;
    legendInfo{li_ii} = 'Hit';
end
hold on
p2=polarplot(pi*meanVectorAngle(graphMiss==1)/180,percent_lick(graphMiss==1),' vr','MarkerFaceColor','r', 'MarkerSize',3);
if ~isempty(p2)
    li_ii=li_ii+1;
    legendInfo{li_ii} = 'Miss';
end
p3=polarplot(pi*meanVectorAngle(graphFA==1)/180,percent_lick(graphFA==1),' vb','MarkerFaceColor','b', 'MarkerSize',3);
if ~isempty(p3)
    li_ii=li_ii+1;
    legendInfo{li_ii} = 'FA';
end
p4=polarplot(pi*meanVectorAngle(graphCR==1)/180,percent_lick(graphCR==1),' vb', 'MarkerFaceColor','b', 'MarkerSize',7);
if ~isempty(p4)
    li_ii=li_ii+1;
    legendInfo{li_ii} = 'CR';
end
legend(legendInfo)
title('Pecent lick as a function of mean phase')

try
    close 7
catch
end

hFig7 = figure(7);
set(hFig6, 'units','normalized','position',[.49 .05 .35 .35])

plot(mod_indx(graphHit==1),percent_lick(graphHit==1),' vr', 'MarkerFaceColor','r', 'MarkerSize',7)
hold on
plot(mod_indx(graphMiss==1),percent_lick(graphMiss==1),' vr','MarkerFaceColor','r', 'MarkerSize',3)
plot(mod_indx(graphFA==1),percent_lick(graphFA==1),' vb','MarkerFaceColor','b', 'MarkerSize',3)
plot(mod_indx(graphCR==1),percent_lick(graphCR==1),' vb', 'MarkerFaceColor','b', 'MarkerSize',7)
legend(legendInfo)
xlabel('Modulation index')
ylabel('Percent lick')
title('Modulation index vs percent lick')

pffft=1




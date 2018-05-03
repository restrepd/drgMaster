function [rho,which_event,no_trials,angleLFPexp,filtLFPexp,filtLFPref]=drgComparePhases(handles)
tic
%This function compares the Hibert transform phases and
%computes the correlation between the filtered LFPs
anglerefLFP = [];
angleLFPexp = [];
filtLFPref=[];
filtLFPexp=[];
delta_phase_timecourse=[];

%Generates a trial per trial phase histogram
sessionNo=handles.sessionNo;
Fs=floor(handles.drg.session(sessionNo).draq_p.ActualRate);
lowF1=handles.peakLowF;
lowF2=handles.peakHighF;
highF1=handles.burstLowF;
highF2=handles.burstHighF;
pad_time=handles.time_pad;
pad_ii=floor(pad_time*Fs);
delta_ii=floor(Fs*(handles.time_end-handles.time_start-2*pad_time));
n_phase_bins=handles.n_phase_bins;

%Enter trials
firstTr=handles.trialNo;
lastTr=handles.lastTrialNo;


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

for trNo=firstTr:lastTr
    
    if handles.displayData==1
        trial_no=trNo
    end
    
    evNo = drgFindEvNo(handles,trNo,sessionNo);
    if evNo~=-1
        excludeTrial=drgExcludeTrialLFP(handles.drg,handles.peakLFPNo,handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);
        
        if excludeTrial==0
            
            %Note: handles.peakLFPNo is the reference LFP
            
            [LFPref, trialNo, can_read1] = drgGetTrialLFPData(handles, handles.peakLFPNo, evNo, handles.evTypeNo, handles.time_start, handles.time_end);
            [LFPexp, trialNo, can_read2] = drgGetTrialLFPData(handles, handles.burstLFPNo, evNo, handles.evTypeNo, handles.time_start, handles.time_end);
            
            if (can_read1==1)&(can_read2==1)
                no_trials=no_trials+1;
                time(no_trials)=handles.drg.session(sessionNo).trial_start(trialNo);
                which_trial(no_trials)=1;
                perCorr_per_histo(no_trials)=50;
                
                %Get theta phase
                bpFiltLFPref = designfilt('bandpassiir','FilterOrder',20, ...
                    'HalfPowerFrequency1',lowF1,'HalfPowerFrequency2',lowF2, ...
                    'SampleRate',Fs);
                thisfiltLFPref=filtfilt(bpFiltLFPref,LFPref);
                thisanglerefLFP = angle(hilbert(thisfiltLFPref)); % phase modulation of theta amplitude
                dectha=decimate(thisanglerefLFP(pad_ii+1:pad_ii+delta_ii),20);
                anglerefLFP(no_trials,1:length(dectha)) = dectha;
                decthl=decimate(thisfiltLFPref(pad_ii+1:pad_ii+delta_ii),20);
                filtLFPref(no_trials,1:length(decthl))=decthl;
                
                %Get LFPexp phase
                bpFiltLFPexp = designfilt('bandpassiir','FilterOrder',20, ...
                    'HalfPowerFrequency1',highF1,'HalfPowerFrequency2',highF2, ...
                    'SampleRate',Fs);
                thisfiltLFPexp=filtfilt(bpFiltLFPexp,LFPexp);
                thisangleLFPexp = angle(hilbert(thisfiltLFPexp)); % LFPexp phase
             
                decale=decimate(thisangleLFPexp(pad_ii+1:pad_ii+delta_ii),20);
                angleLFPexp(no_trials,1:length(decale)) = decale;
                decfle=decimate(thisfiltLFPexp(pad_ii+1:pad_ii+delta_ii),20);
                filtLFPexp(no_trials,1:length(decfle))=decfle;
                
                
                rho(no_trials)=corr(thisfiltLFPref(pad_ii+1:pad_ii+delta_ii)',thisfiltLFPexp(pad_ii+1:pad_ii+delta_ii)');
                
                %Delta phase
                delta_phase_timecourse(no_trials,1:length(decthl))=decimate(thisangleLFPexp(pad_ii+1:pad_ii+delta_ii)-thisanglerefLFP(pad_ii+1:pad_ii+delta_ii),20);
                
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
    end
    
    %end %if eventstamps...
end %for evNo

delta_t=0.2;
edges=[-2*pi():4*pi()/50:2*pi()];
delta_phase_t_hist=zeros(floor((handles.time_end-handles.time_start-2*pad_time)/delta_t),50);
delta_ii=floor(delta_t*Fs/20);
time=handles.time_start+pad_time;

for jj=1:floor((handles.time_end-handles.time_start-2*pad_time)/delta_t)
    t(jj)=handles.time_start+pad_time+(delta_t/2)+(jj-1)*delta_t;
    for trNo=1:no_trials
        these_phases=delta_phase_timecourse(trNo,(jj-1)*delta_ii+1:jj*delta_ii);
        delta_phase_t_hist(jj,1:50)=delta_phase_t_hist(jj,1:50)+histcounts(these_phases,edges)/(delta_ii*no_trials);
    end
end

toc

if handles.displayData==1
    try
        close 2
    catch
    end
    
    hFig2 = figure(2);
    set(hFig2, 'units','normalized','position',[.75 .1 .25 .3])
    
    edges=[-2*pi():4*pi()/50:2*pi()];
    histogram(angleLFPexp(:)-anglerefLFP(:),edges,'Normalization','probability')
    title('Histogram of the difference in phase between LFPexp and LFPref')
    ylabel('Probability')
    xlabel('Delta phase LFPexp - LFPref (radians)')
    
    try
        close 3
    catch
    end
    
    hFig3 = figure(3);
    set(hFig3, 'units','normalized','position',[.01 .1 .67 .3])
    
    
   
    
    
    phase=[-2*pi()+(4*pi()/100):4*pi()/50:2*pi()-(4*pi()/100)]';
    drg_pcolor(repmat(t,50,1),repmat(phase,1,length(t)),delta_phase_t_hist')
    colormap jet
    shading flat
    
    min_prob=prctile(delta_phase_t_hist(:),1);
    max_prob=prctile(delta_phase_t_hist(:),99);
    % caxis([min_prob    max_prob])
    caxis([0 0.2])
    
    xlabel('Time (sec)')
    ylabel('Radians');
    title('Timecourse for histogram of delta phase LFPexp - LFPref (radians)')
    
    try
        close 4
    catch
    end
    
    hFig4 = figure(4);
    set(hFig4, 'units','normalized','position',[.69 .1 .05 .3])
    
    prain=[min_prob:(max_prob-min_prob)/99:max_prob];
    pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
    colormap jet
    shading interp
    ax=gca;
    set(ax,'XTickLabel','')
    
    try
        close 1
    catch
    end
    
    hFig1 = figure(1);
    set(hFig1, 'units','normalized','position',[.01 .55 .25 .3])
    
    plot(anglerefLFP(:),angleLFPexp(:),'ob','MarkerSize',0.6)
    xlabel(['Phase LFPref electrode #' num2str(handles.peakLFPNo)])
    ylabel(['Phase LFPexp electrode #' num2str(handles.burstLFPNo)])
    xlim([-pi() pi()])
    ylim([-pi() pi()])
    
    try
        close 5
    catch
    end
    
    hFig5 = figure(5);
    set(hFig5, 'units','normalized','position',[.26 .55 .25 .3])
    
    plot(filtLFPref(:),filtLFPexp(:),'ob','MarkerSize',0.6)
    xlabel(['LFPref electrode #' num2str(handles.peakLFPNo)])
    ylabel(['LFPexp electrode #' num2str(handles.burstLFPNo)])
    pct99=prctile([filtLFPref(:); filtLFPexp(:)],99);
    pct1=prctile([filtLFPref(:); filtLFPexp(:)],1);
    xlim([pct1 pct99])
    ylim([pct1 pct99])
    
    try
        close 6
    catch
    end
    
    hFig6 = figure(6);
    set(hFig6, 'units','normalized','position',[.01 .55 .25 .3])
    
    delta_t=0.2;
    delta_ii=floor(delta_t*(Fs/20));
    rho_timecourse=zeros(1,floor((handles.time_end-handles.time_start-2*pad_time)/delta_t));
    
    
    for jj=1:floor((handles.time_end-handles.time_start-2*pad_time)/delta_t)
        t(jj)=handles.time_start+pad_time+(delta_t/2)+(jj-1)*delta_t;
        these_filtLFPexp=[];
        these_filtLFPref=[];
        for trNo=1:no_trials
            these_filtLFPexp=[these_filtLFPexp filtLFPexp(trNo,(jj-1)*delta_ii+1:jj*delta_ii)];
            these_filtLFPref=[these_filtLFPref filtLFPref(trNo,(jj-1)*delta_ii+1:jj*delta_ii)];
        end
        
        rho_timecourse(jj)=corr(these_filtLFPref',these_filtLFPexp');
        
    end
    
    plot(t,rho_timecourse,'-o')
    
    xlabel('Time (sec)')
    ylabel('Rho');
    title('Timecourse for correlation coefficient LFPexp vs LFPref')
    ylim([0 1])
    
    fprintf(1, 'Mean correlation coefficient=%d\n',mean(rho));
    
    pffft=1;
end







function drgdFvsBFREncRetr(handles)
%PSTH plot the encoding and retrieval segments

startBasal=-2.5;
endBasal=0;
midBasal=-1.25;
startOdor=0;
endOdor=2.5;

sessionNo=handles.drg.unit(handles.unitNo).sessionNo;


[perCorr, encoding_trials, retrieval_trials, encoding_this_evTypeNo,retrieval_this_evTypeNo]=drgFindEncRetr(handles);
%Enter unit and event
%unitNo=2;

try
    close 1
catch
end

hFig1=figure(1);
set(hFig1, 'units','normalized','position',[.7 .1 .3 .3])

%Plot the percent correct
trials=1:length(perCorr);
plot(trials,perCorr,'o','MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7])
hold on
plot(trials(encoding_trials),perCorr(encoding_trials),'ob')
plot(trials(retrieval_trials),perCorr(retrieval_trials),'or')
ylim([40 110]);
xlabel('Trial #')
ylabel('Percent correct')
title('% correct vs trial number, blue=encoding, red=retrieval')

%Plot BFR pre1 vs pre 2 for all trials in the session
noTrials=0;
preBFR1=[];
preBFR2=[];
odorOn=2;
spike_times=handles.drg.unit(handles.unitNo).spike_times;
for evNo=1:handles.drg.session(sessionNo).events(odorOn).noTimes
    
    %evNo
    excludeTrial=drgExcludeTrial(handles.drg,handles.drg.unit(handles.unitNo).channel,handles.drg.session(sessionNo).events(odorOn).times(evNo),sessionNo);
    
    if excludeTrial==0
        
        noTrials=noTrials+1;
        
        preBFR1spikes=[];
        preBFR1spikes=(spike_times>handles.drg.session(sessionNo).events(odorOn).times(evNo)+startBasal)&...
            (spike_times<=handles.drg.session(sessionNo).events(odorOn).times(evNo)+midBasal);
        preBFR1(noTrials)=sum(preBFR1spikes)/(midBasal-startBasal);
        
        preBFR2spikes=[];
        preBFR2spikes=(spike_times>handles.drg.session(sessionNo).events(odorOn).times(evNo)+midBasal)&...
            (spike_times<=handles.drg.session(sessionNo).events(odorOn).times(evNo)+endBasal);
        preBFR2(noTrials)=sum(preBFR2spikes)/(endBasal-midBasal);
        
    end
    %end
    %end %if eventstamps...
end %for evNo

try
    close(2)
catch
end
hFig2=figure(2);
set(hFig2, 'units','normalized','position',[.37 .1 .3 .3])
hold on


minBFR1=min(preBFR1);
maxBFR1=max(preBFR1);
minBFR2=min(preBFR2);
maxBFR2=max(preBFR2);
xlim([minBFR1-0.2*(maxBFR1-minBFR1) maxBFR1+0.2*(maxBFR1-minBFR1)])
ylim([minBFR2-0.2*(maxBFR2-minBFR2) maxBFR2+0.2*(maxBFR2-minBFR2)])

beta_BFR12=nlinfit(preBFR1',preBFR2',@dr_fitline,[0 1]);


[rho pval_rho]=corr(preBFR1',preBFR2');

%slope_12(unitPairNo)=beta_BFR12(1);
fitBFR2 = dr_fitline(beta_BFR12,[minBFR1-0.2*(maxBFR1-minBFR1) maxBFR1+0.2*(maxBFR1-minBFR1)]);
plot([minBFR1-0.2*(maxBFR1-minBFR1) maxBFR1+0.2*(maxBFR1-minBFR1)],fitBFR2,'-b')

plot(preBFR1,preBFR2,'ob')
%plot(preBFR1,postBFR2,'or')
title('BFR2 vs BFR1')
xlabel('BFR1, Hz')
ylabel('BFR2, Hz')


if pval_rho<0.05
    text(minBFR1-0.05*(maxBFR1-minBFR1),minBFR2-0.05*(maxBFR2-minBFR2),['rho= ' num2str(rho)],'Color','r')
else
    text(minBFR1-0.05*(maxBFR1-minBFR1),minBFR2-0.05*(maxBFR2-minBFR2),['rho= ' num2str(rho)],'Color','b')
end


%Do encoding
bin_size=0.1;
nobins=floor((endOdor-startBasal)/bin_size);
startBasalBin=1;
endBasalBin=24;
startOdorBin=25;
endOdorBin=50;
itime=1:nobins;
itime=itime+fix(handles.drg.time_pre/bin_size);
time=double(itime)*bin_size;
noTrials_encoding=0;
spike_times=[];
spike_times=handles.drg.unit(handles.unitNo).spike_times;
PSTHencoding=[];


for evNo=encoding_this_evTypeNo
    
    %evNo
    excludeTrial=drgExcludeTrial(handles.drg,handles.drg.unit(handles.unitNo).channel,handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);
    
    if excludeTrial==0
        
        noTrials_encoding=noTrials_encoding+1;
        these_spikes=(spike_times>handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo)+startBasal)&...
            (spike_times<handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo)+endOdor);
        these_spike_times=spike_times(these_spikes)-handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo);
        PSTHencoding(noTrials_encoding,1:nobins)=0;
        for spk=1:length(these_spike_times)
            this_bin=ceil((these_spike_times(spk)-startBasal)/bin_size);
            PSTHencoding(noTrials_encoding,this_bin)=PSTHencoding(noTrials_encoding,this_bin)+1;
        end %for spk
        %Calculate the p values per trial
        [h p_val_within_pre_vs_post_enc(noTrials_encoding)]=ttest2(PSTHencoding(noTrials_encoding,startBasalBin:endBasalBin),PSTHencoding(noTrials_encoding,startOdorBin:endOdorBin));
    end
    %end
    %end %if eventstamps...
end %for evNo

p_val_enc_FDR=drsFDRpval(p_val_within_pre_vs_post_enc);

dF_enc=(mean(PSTHencoding(:,startOdorBin:endOdorBin),2)-mean(PSTHencoding(:,startBasalBin:endBasalBin),2))/bin_size;
BFRenc=mean(PSTHencoding(:,startBasalBin:endBasalBin),2)/bin_size;

minBFRenc=min(BFRenc);
maxBFRenc=max(BFRenc);
maxdFenc=max(dF_enc);
mindFenc=min(dF_enc);

beta_enc=nlinfit(BFRenc,dF_enc,@dr_fitline,[0 1]);


%Now plot the dF vs BFR
try
    close(3)
catch
end
hFig3=figure(3);
set(hFig3, 'units','normalized','position',[.05 .6 .3 .3])

plot(BFRenc(p_val_within_pre_vs_post_enc<=p_val_enc_FDR),dF_enc(p_val_within_pre_vs_post_enc<=p_val_enc_FDR),'o','MarkerEdgeColor','b','MarkerFaceColor','b')
hold on
plot(BFRenc(p_val_within_pre_vs_post_enc>p_val_enc_FDR),dF_enc(p_val_within_pre_vs_post_enc>p_val_enc_FDR),'o','MarkerEdgeColor','b','MarkerFaceColor','none')
plot([minBFRenc-0.2*(maxBFRenc-minBFRenc) maxBFRenc+0.2*(maxBFRenc-minBFRenc)],[0 0],'-k')
plot([minBFRenc-0.2*(maxBFRenc-minBFRenc) maxBFRenc+0.2*(maxBFRenc-minBFRenc)],dr_fitline(beta_enc,[minBFRenc-0.2*(maxBFRenc-minBFRenc) maxBFRenc+0.2*(maxBFRenc-minBFRenc)]),'-b')
plot([mean(BFRenc) mean(BFRenc)],[mindFenc-0.2*(maxdFenc-mindFenc) maxdFenc+0.2*(maxdFenc-mindFenc)],'-r')

ylim([mindFenc-0.2*(maxdFenc-mindFenc) maxdFenc+0.2*(maxdFenc-mindFenc)])
xlim([minBFRenc-0.2*(maxBFRenc-minBFRenc) maxBFRenc+0.2*(maxBFRenc-minBFRenc)])

title(['deltaFR vs. BFR, encoding ' handles.drg.session.eventlabels{handles.evTypeNo}])
ylabel('delta FR(Hz)')
xlabel('BFR (Hz)')

%Do retrieval
noTrials_retreival=0;
PSTHretreival=[];


for evNo=retrieval_this_evTypeNo
    
    %evNo
    excludeTrial=drgExcludeTrial(handles.drg,handles.drg.unit(handles.unitNo).channel,handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);
    
    if excludeTrial==0
        
        noTrials_retreival=noTrials_retreival+1;
        these_spikes=(spike_times>handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo)+startBasal)&...
            (spike_times<handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo)+endOdor);
        these_spike_times=spike_times(these_spikes)-handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo);
        PSTHretreival(noTrials_retreival,1:nobins)=0;
        for spk=1:length(these_spike_times)
            this_bin=ceil((these_spike_times(spk)-startBasal)/bin_size);
            PSTHretreival(noTrials_retreival,this_bin)=PSTHretreival(noTrials_retreival,this_bin)+1;
        end %for spk
        %Calculate the p values per trial
        [h p_val_within_pre_vs_post_retr(noTrials_retreival)]=ttest2(PSTHretreival(noTrials_retreival,startBasalBin:endBasalBin),PSTHretreival(noTrials_retreival,startOdorBin:endOdorBin));
    end
    %end
    %end %if eventstamps...
end %for evNo


p_val_retr_FDR=drsFDRpval(p_val_within_pre_vs_post_retr);

dF_retr=(mean(PSTHretreival(:,startOdorBin:endOdorBin),2)-mean(PSTHretreival(:,startBasalBin:endBasalBin),2))/bin_size;
BFRretr=mean(PSTHretreival(:,startBasalBin:endBasalBin),2)/bin_size;

minBFRretr=min(BFRretr);
maxBFRretr=max(BFRretr);
maxdFretr=max(dF_retr);
mindFretr=min(dF_retr);

beta_retr=nlinfit(BFRretr,dF_retr,@dr_fitline,[0 1]);


%Now plot the dF vs BFR
try
    close(4)
catch
end
hFig4=figure(4);
set(hFig4, 'units','normalized','position',[.05 .1 .3 .3])

plot(BFRretr(p_val_within_pre_vs_post_retr<=p_val_retr_FDR),dF_retr(p_val_within_pre_vs_post_retr<=p_val_retr_FDR),'o','MarkerEdgeColor','b','MarkerFaceColor','b')
hold on
plot(BFRretr(p_val_within_pre_vs_post_retr>p_val_retr_FDR),dF_retr(p_val_within_pre_vs_post_retr>p_val_retr_FDR),'o','MarkerEdgeColor','b','MarkerFaceColor','none')
plot([minBFRretr-0.2*(maxBFRretr-minBFRretr) maxBFRretr+0.2*(maxBFRretr-minBFRretr)],[0 0],'-k')
plot([minBFRretr-0.2*(maxBFRretr-minBFRretr) maxBFRretr+0.2*(maxBFRretr-minBFRretr)],dr_fitline(beta_retr,[minBFRretr-0.2*(maxBFRretr-minBFRretr) maxBFRretr+0.2*(maxBFRretr-minBFRretr)]),'-b')
plot([mean(BFRretr) mean(BFRretr)],[mindFretr-0.2*(maxdFretr-mindFretr) maxdFretr+0.2*(maxdFretr-mindFretr)],'-r')

ylim([mindFretr-0.2*(maxdFretr-mindFretr) maxdFretr+0.2*(maxdFretr-mindFretr)])
xlim([minBFRretr-0.2*(maxBFRretr-minBFRretr) maxBFRretr+0.2*(maxBFRretr-minBFRretr)])

title(['deltaFR vs. BFR, retreival ' handles.drg.session.eventlabels{handles.evTypeNo}])
ylabel('delta FR(Hz)')
xlabel('BFR (Hz)')

%Plot the PSTH for each trial
try
    close(8)
catch
end
hFig8=figure(8);
set(hFig8, 'units','normalized','position',[.4 .05 .25 .7])
hold on

to_sort_enc=[[1:length(BFRenc)]' BFRenc];
sorted_enc=sortrows(to_sort_enc,2);

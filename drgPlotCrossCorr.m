function handles=drgPlotCrossCorr(handles)
%Plot autocorrelaogram

noShuffles=20;
sessionNo=handles.drg.unit(handles.unitNo).sessionNo;

%Enter trials
firstTr=handles.trialNo;
lastTr=handles.lastTrialNo;


auto_width=handles.corr_window;
nobins=21;
bin_size=2*handles.corr_window/nobins;
delta_times=[-auto_width+bin_size/2:bin_size:auto_width-bin_size/2];
noTrials=0;
spike_times=[];
spike_times=handles.drg.unit(handles.unitNo).spike_times;
spike_times2=handles.drg.unit(handles.unitNo2).spike_times;
Auto=zeros(1,nobins);
no_comp_spikes=0;
shAuto=zeros(1,nobins);
no_comp_spikes_sh=0;
time_start=handles.time_start+handles.time_pad;
time_end=handles.time_end-handles.time_pad;
synch_spikes=[];
synch_spikes_sh=[];
ref_spikes=0;


for trNo=firstTr:lastTr
    
    if handles.save_drgb==0
        trial_no=trNo
    end
    
    evNo = drgFindEvNo(handles,trNo,sessionNo);
    
    if evNo~=-1
        
        %evNo
        excludeTrial=drgExcludeTrial(handles.drg,handles.drg.unit(handles.unitNo).channel,handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);
        
        if excludeTrial==0
            
            noTrials=noTrials+1;
            synch_spikes(noTrials)=0;
            synch_spikes_sh(noTrials)=0;
            these_spikes=(spike_times>handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo)+time_start+auto_width)&...
                (spike_times<=handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo)+time_end-auto_width);
            these_spike_times=spike_times(these_spikes)-(handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo)+time_start);
            ref_spikes=ref_spikes+sum(these_spikes);
            these_spikes2=(spike_times2>handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo)+time_start+auto_width)&...
                (spike_times2<=handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo)+time_end-auto_width);
            these_spike_times2=spike_times2(these_spikes2)-(handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo)+time_start);
            
            for spkref=1:length(these_spike_times)
                for spk=1:length(these_spike_times2)
                    
                    deltat=these_spike_times2(spk)-these_spike_times(spkref);
                    if abs(deltat)<auto_width+bin_size/2
                        %The interspike interval is within auto_width
                        synch_spikes(noTrials)=synch_spikes(noTrials)+1;
                        this_bin=find((delta_times-(bin_size/2)<=deltat)&(delta_times+(bin_size/2)>deltat));
                        if (this_bin>0)&(this_bin<=nobins) 
                            Auto(1,this_bin)=Auto(1,this_bin)+1;
                        end
                    end
                    
                end
            end %for spkref
            no_comp_spikes=no_comp_spikes+length(these_spike_times2)-1;
            
            %Now do random shuffled spikes
            for noS=1:noShuffles
                shuf_spike_times=(time_end-time_start-2*+auto_width)*rand(1,length(these_spike_times))+auto_width;
                shuf_spike_times2=(time_end-time_start-2*+auto_width)*rand(1,length(these_spike_times2))+auto_width;
                for spkref=1:length(these_spike_times)
                    for spk=1:length(these_spike_times2)
                        
                        deltat=shuf_spike_times2(spk)-shuf_spike_times(spkref);
                        if abs(deltat)<auto_width
                            synch_spikes_sh(noTrials)=synch_spikes_sh(noTrials)+1;
                            this_bin=find((delta_times-(bin_size/2)<=deltat)&(delta_times+(bin_size/2)>deltat));
                            if (this_bin>0)&(this_bin<=nobins)
                                shAuto(1,this_bin)=shAuto(1,this_bin)+1;
                            end
                            
                            
                        end
                        
                    end
                end %for spkref
                no_comp_spikes_sh=no_comp_spikes_sh+length(these_spike_times2)-1;
            end
        end
        %end
        %end %if eventstamps...
    end %if evNo
end %for trNo=

no_trials_included=noTrials;

if no_comp_spikes~=0
    Auto=Auto/no_comp_spikes;
end

if no_comp_spikes_sh~=0
    shAuto=shAuto/no_comp_spikes_sh;
    synch_spikes_sh=synch_spikes_sh/no_comp_spikes_sh;
end

if handles.displayData==1
    %Now plot the crosscorrelogram
    try
        close 1
    catch
    end
    
    %Plot the timecourse
    hFig1 = figure(1);
    set(hFig1, 'units','normalized','position',[.02 .4 .5 .3])
    bar(delta_times,Auto,'b');

    
    %Now plot the shuffled autocorrelogram - random
    hold on
    bar(delta_times,shAuto,'FaceColor','none','EdgeColor','r');
    title(['Cross Correlogram for ' handles.drg.session.eventlabels{handles.evTypeNo}])
    ylabel('Correlation coefficient')
    xlabel('delta time (sec)')

    y_max=max(Auto);
    if y_max~=0
        ylim([0 1.2*y_max])
    else
        ylim([0 1.2])
    end
end

%Calculate the modulation index defined by Tort et al J Neurophysiol 104: 1195?1210, 2010
%Note that the pvalue for Tort et al is the same as phase_histo
mean_prob=mean(Auto)*ones(1,length(Auto));
Auto_non_zero=Auto+0.00000000001;
DKL=sum(Auto_non_zero(1:end).*log(Auto_non_zero(1:end)./mean_prob(1:end)));
MI_Tort=DKL/log(length(Auto));
skew=skewness(Auto);
p_value=ranksum(synch_spikes_sh,synch_spikes);
percent_synch=100*(sum(synch_spikes)-sum(synch_spikes_sh))/ref_spikes;

%Save the data if requested by drgb
if handles.save_drgb==1
    handles.drgb.corrpair(handles.drgb.corrpair_no).delta_times=delta_times;
    handles.drgb.corrpair(handles.drgb.corrpair_no).correlogram=Auto;
    handles.drgb.corrpair(handles.drgb.corrpair_no).MI_Tort=MI_Tort;
    handles.drgb.corrpair(handles.drgb.corrpair_no).skewness=skew;
    handles.drgb.corrpair(handles.drgb.corrpair_no).p_value=p_value;
    handles.drgb.corrpair(handles.drgb.corrpair_no).percent_synch=percent_synch;
    
    %Reference unit
    handles.drgb.corrpair(handles.drgb.corrpair_no).ref_unit=handles.unitNo;
    handles.drgb.corrpair(handles.drgb.corrpair_no).ref_u_single=handles.drg.unit(handles.unitNo).SingleUnit;
    handles.drgb.corrpair(handles.drgb.corrpair_no).ref_u_perViol=handles.drg.unit(handles.unitNo).perViol;
    handles.drgb.corrpair(handles.drgb.corrpair_no).ref_u_Lratio=handles.drg.unit(handles.unitNo).Lratio;
    
    %Partner unit
    handles.drgb.corrpair(handles.drgb.corrpair_no).partner_unit=handles.unitNo2;
    handles.drgb.corrpair(handles.drgb.corrpair_no).partner_u_single=handles.drg.unit(handles.unitNo2).SingleUnit;
    handles.drgb.corrpair(handles.drgb.corrpair_no).partner_u_perViol=handles.drg.unit(handles.unitNo2).perViol;
    handles.drgb.corrpair(handles.drgb.corrpair_no).partner_u_Lratio=handles.drg.unit(handles.unitNo2).Lratio;
    
end
 
pffft=1;

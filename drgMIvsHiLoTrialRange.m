function drgMIvsHiLoTrialRange(handles)

deltaLowF_PAC=handles.deltaLowF_PAC;
deltaHighF_PAC=handles.deltaHighF_PAC;
bandwidth_lowF=handles.bandwidth_lowF;
bandwidth_highF=handles.bandwidth_highF;

%Generates a comodulogram for the encoding and retrieval segments
sessionNo=handles.sessionNo;
Fs=handles.drg.session(sessionNo).draq_p.ActualRate;
lowF1=handles.peakLowF-(handles.bandwidth_lowF/2);
lowF2=handles.peakHighF-(handles.bandwidth_lowF/2);
highF1=handles.burstLowF-(handles.bandwidth_highF/2);
highF2=handles.burstHighF-(handles.bandwidth_highF/2);
pad_time=handles.time_pad;
n_phase_bins=handles.n_phase_bins;

%Enter trials
firstTr=handles.trialNo;
lastTr=handles.lastTrialNo;




%Now get the phase of gamma bursts witin theta
no_encoding_trials=0;
no_retrieval_trials=0;
comodulogram_enc=[];
compdulogram_retr=[];
phase_histo_enc=zeros(length([lowF1:deltaLowF_PAC:lowF2]),length(highF1:deltaHighF_PAC:highF2),51);
phase_histo_retr=zeros(length([lowF1:deltaLowF_PAC:lowF2]),length(highF1:deltaHighF_PAC:highF2),51);

lowF6ii=[];
highF65ii=[];

for trNo=firstTr:lastTr
    
    if handles.save_drgb==0
        trialNo=trNo
    end
    
    evNo = drgFindEvNo(handles,trNo,sessionNo);
    
    if evNo~=-1
        
        excludeTrial=drgExcludeTrialLFP(handles.drg,handles.peakLFPNo,handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);
        
        if excludeTrial==0
            
            
            [LFPlow, trialNo, can_read1] = drgGetTrialLFPData(handles, handles.peakLFPNo, evNo, handles.evTypeNo, handles.time_start, handles.time_end);
            [LFPhigh, trialNo, can_read2] = drgGetTrialLFPData(handles, handles.burstLFPNo, evNo, handles.evTypeNo, handles.time_start, handles.time_end);
            
            if (can_read1==1)&(can_read2==1)
                
                lowFii=0;
                for lowF=lowF1:deltaLowF_PAC:lowF2
                    lowFii=lowFii+1;
                    highFii=0;
                    for highF=highF1:deltaHighF_PAC:highF2
                        
                        [meanVectorLength, meanVectorAngle, peakAngle, mod_indx, phase, phase_histo, theta_wave]=drgGetThetaAmpPhase(LFPlow, LFPhigh,Fs,lowF,lowF+bandwidth_lowF,highF,highF+bandwidth_highF,pad_time,n_phase_bins,handles.which_method);
                        
                        highFii=highFii+1;
                        
                        %                         if (lowF==6)&(highF==65)
                        %                             lowF6ii=[lowF6ii lowFii];
                        %                             highF65ii=[highF65ii highFii];
                        %                         end
                        %
                        
                        if (lowF==lowF1)&(highF==highF1)
                            no_encoding_trials=no_encoding_trials+1;
                        end
                        comodulogram_enc(no_encoding_trials,lowFii,highFii)=mod_indx;
                        add_hist=zeros(length([lowF1:deltaLowF_PAC:lowF2]),length(highF1:deltaHighF_PAC:highF2),51);
                        add_hist(lowFii,highFii,:)=phase_histo;
                        phase_histo_enc=phase_histo_enc+add_hist;
                        
                        
                        
                        
                        
                    end
                    
                end
                
                
            end
        end
    end
    %end
    %end %if eventstamps...
end %for evNo



%IMPORTANT: pcolor does not plot the last column/row
%Because of this we have to add one row/column
lowFii=0;
for lowF=lowF1:deltaLowF_PAC:lowF2+deltaLowF_PAC
    lowFii=lowFii+1;
    highFii=0;
    for highF=highF1:deltaHighF_PAC:highF2+deltaHighF_PAC
        highFii=highFii+1;
        lowF_matrix(lowFii,highFii)=lowF+(bandwidth_lowF/2)-(deltaLowF_PAC/2);
        highF_matrix(lowFii,highFii)=highF+(bandwidth_highF/2);
    end
end


szcomenc=size(comodulogram_enc);


mean_com_enc=zeros(szcomenc(2)+1,szcomenc(3)+1);
mean_com_enc(1:szcomenc(2),1:szcomenc(3))=mean(comodulogram_enc,1);



max_com=max(max(mean_com_enc(1:szcomenc(2),1:szcomenc(3))));
min_com=min(min(mean_com_enc(1:szcomenc(2),1:szcomenc(3))));

try
    close 1
catch
end

hFig1 = figure(1);
set(hFig1, 'units','normalized','position',[.01 .1 .23 .8])

pcolor(lowF_matrix,highF_matrix,mean_com_enc)
colormap jet
caxis([min_com max_com]);
shading interp
xlabel('Frequency for phase (Hz)')
ylabel('Frequency for amplitude (Hz)');
title(['Phase-amplitude comodulogram ' handles.drg.session(1).draq_d.eventlabels{handles.evTypeNo}])


try
    close 2
catch
end

hFig2 = figure(2);
set(hFig2, 'units','normalized','position',[.49 .1 .05 .3])

prain=[min_com:(max_com-min_com)/99:max_com];
pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
colormap jet
shading interp
ax=gca;
set(ax,'XTickLabel','')

%Calculate MI for the average phase histogram
%Calculate the modulation index defined by Tort et al J Neurophysiol 104: 1195?1210, 2010
%Note that the pvalue for Tort et al is the same as phase_histo
phase_histo_enc=phase_histo_enc/no_encoding_trials;
szphe=size(phase_histo_enc);
MI_Tort_enc=zeros(lowFii,highFii);
for lowFkk=1:lowFii
    for highFkk=1:highFii
        if (lowFkk~=lowFii)&(highFkk~=highFii)
            mean_prob_enc=zeros(1,szphe(3));
            mean_prob_enc(1,:)=ones(1,szphe(3))*mean(phase_histo_enc(lowFkk,highFkk,:));
            phe=zeros(1,szphe(3)-1);
            phe(1,:)=phase_histo_enc(lowFkk,highFkk,1:end-1);
            DKL_enc=sum(phe.*log(phe./mean_prob_enc(1,1:end-1)));
            MI_Tort_enc(lowFkk,highFkk)=DKL_enc/log(n_phase_bins);
        else
            MI_Tort_enc(lowFkk,highFkk)=0;
        end
    end
end



max_MI=max(max(MI_Tort_enc(1:szcomenc(2),1:szcomenc(3))));
min_MI=min(min(mean_com_enc(1:szcomenc(2),1:szcomenc(3))));


try
    close 3
catch
end

hFig3 = figure(3);
set(hFig3, 'units','normalized','position',[.01 .05 .23 .8])

pcolor(lowF_matrix,highF_matrix,MI_Tort_enc)
colormap jet
caxis([min_MI max_MI]);
shading interp
xlabel('Frequency for phase (Hz)')
ylabel('Frequency for amplitude (Hz)');
title(['Phase-amplitude MI comodulogram/mean phase' handles.drg.session(1).draq_d.eventlabels{handles.evTypeNo}])



try
    close 4
catch
end

hFig4 = figure(4);
set(hFig4, 'units','normalized','position',[.49 .05 .05 .3])

prain=[min_MI:(max_MI-min_MI)/99:max_MI];
pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
colormap jet
shading interp
ax=gca;
set(ax,'XTickLabel','')

pffft=1






function drgLFPimgcohDeltaallPairs(handles)
%Generates a timecourse of the imaginary coherence between two LFP channels
 
current_peakLFPNo=handles.peakLFPNo;
current_burstLFPNo=handles.burstLFPNo;

deltaCxy=[];
iiCxy=0;
t_pre=[-2 0];
t_post=[0.5 2.5];

for this_peakLFPNo=1:8
    for this_burstLFPNo=9:16
        handles.peakLFPNo=this_peakLFPNo;
        handles.burstLFPNo=this_burstLFPNo;
        [t,f, all_Cxy_timecourse, this_trialNo]=drgGetLFPimgCoherenceForThisEvTypeNo(handles);
        
        pre_t=(t>=t_pre(1))&(t<=t_pre(2));
        Cxy_pre_per_trial=zeros(1,size(all_Cxy_timecourse,1));
        Cxy_pre_per_trial(:,:)=mean(mean(all_Cxy_timecourse(:,(f>=6)&(f<=14),pre_t),2),3);
        
        post_t=(t>=t_post(1))&(t<=t_post(2));
        Cxy_post_per_trial=zeros(1,size(all_Cxy_timecourse,1));
        Cxy_post_per_trial(:,:)=mean(mean(all_Cxy_timecourse(:,(f>=6)&(f<=14),post_t),2),3);
        
        iiCxy=iiCxy+1;
        deltaCxy(iiCxy)=mean(Cxy_post_per_trial-Cxy_pre_per_trial);
    end
end

handles.peakLFPNo=current_peakLFPNo;
handles.burstLFPNo=current_burstLFPNo;

try
    close(1)
catch
end
hFig1 = figure(1);
set(hFig1, 'units','normalized','position',[.07 .1 .3 .3])


  
edges=[-1:0.05:1];
histogram(deltaCxy,edges)


hold on
plot([mean(deltaCxy) mean(deltaCxy)],[0 12],'-k','LineWidth',2)
ax=gca;ax.LineWidth=3;
ylim([0 14])
xlim([-0.6 0.6])
xlabel('delta Imaginary Coherence')
ylabel('Number')



fprintf(1, ['\n\nMean delta imaginary coherence %d\n'],mean(deltaCxy))


pffft=1

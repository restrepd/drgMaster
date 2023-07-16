function drgLFPimgcohspectTimecourse(handles)
%Generates a timecourse of the imaginary coherence between two LFP channels

[t,f, all_Cxy_timecourse, this_trialNo]=drgGetLFPimgCoherenceForThisEvTypeNo(handles);


freq=f';

Cxy_timecourse(:,:)=mean(all_Cxy_timecourse,1);

if handles.autoscale==1
    maxCxy=prctile(Cxy_timecourse(:),99);
    minCxy=prctile(Cxy_timecourse(:),1);
else
    maxCxy=handles.maxLogP;
    minCxy=handles.minLogP;
end

if minCxy==maxCxy
    minCxy=maxCxy-0.01;
end
figNo=0;

%Plot the timecourse
figNo=figNo+1;
try
    close(figNo)
catch
end
hFig1 = figure(figNo);
set(hFig1, 'units','normalized','position',[.07 .1 .75 .3])


drg_pcolor(repmat(t',length(freq),1)',repmat(f,length(t),1),Cxy_timecourse(1:length(freq),1:length(t))')


colormap fire
shading interp
caxis([minCxy maxCxy]);
xlabel('Time (sec)')
ylabel('Frequency (Hz)');
title(['Coherence timecourse ' handles.drg.session(1).draq_d.eventlabels{handles.evTypeNo}])

figNo=figNo+1;
try
    close(figNo)
catch
end
hFig1 = figure(figNo);
set(hFig1, 'units','normalized','position',[.83 .1 .05 .3])

prain=[minCxy:(maxCxy-minCxy)/99:maxCxy];
drg_pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
colormap fire
shading interp
ax=gca;
set(ax,'XTickLabel','')
% end

%This code is here for Daniels' Figure 1
%Plot the timecourse
figNo=figNo+1;
try
    close(figNo)
catch
end
hFig1 = figure(figNo);
set(hFig1, 'units','normalized','position',[.07 .45 .55 .3])
plot(t',mean(Cxy_timecourse((f>=6)&(f<=14),1:length(t)))','-k','LineWidth',3)
ylim([-1 1])
xlabel('Time (sec)')
ylabel('Coherence')
title('Theta imaginary coherence')

%This code is here for Daniels' Figure 1
%Plot the timecourse
figNo=figNo+1;
try
    close(figNo)
catch
end
hFig1 = figure(figNo);
set(hFig1, 'units','normalized','position',[.07 .45 .55 .3])
plot(t',mean(Cxy_timecourse(:,1:length(t)))','-k','LineWidth',3)
ylim([-1 1])
xlabel('Time (sec)')
ylabel('Coherence')
title('Imaginary coherence for user settings')

%Do a plot of pre-odor vs post-odor if the user chose -2 to 5
t_pre=[-2 0];
t_post=[0.5 2.5];
edges=[-1:0.05:1];
rand_offset=0.8;
ii_rank=0;
input_data=[];

if (sum((t>=t_pre(1))&(t<=t_pre(2)))>5)&(sum((t>=t_post(1))&(t<=t_post(2)))>5)
    figNo=figNo+1;
try
    close(figNo)
catch
end
hFig1 = figure(figNo);
    set(hFig1, 'units','normalized','position',[.07 .7 .15 .3])
    
    ax=gca;ax.LineWidth=3;
    
    hold on
    
    pre_t=(t>=t_pre(1))&(t<=t_pre(2));
    Cxy_pre_per_trial=zeros(1,size(all_Cxy_timecourse,1));
    Cxy_pre_per_trial(:,:)=mean(mean(all_Cxy_timecourse(:,(f>=6)&(f<=14),pre_t),2),3);
    
    bar_offset=1;
    
    bar(bar_offset,mean(Cxy_pre_per_trial),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
    
    
    %Enter the data for t-test/ranksum
    ii_rank=ii_rank+1;
    input_data(ii_rank).data=Cxy_pre_per_trial;
    input_data(ii_rank).description='pre odor';
    
    
    post_t=(t>=t_post(1))&(t<=t_post(2));
    Cxy_post_per_trial=zeros(1,size(all_Cxy_timecourse,1));
    Cxy_post_per_trial(:,:)=mean(mean(all_Cxy_timecourse(:,(f>=6)&(f<=14),post_t),2),3);
    
    bar_offset=bar_offset+1;
    
    bar(bar_offset,mean(Cxy_post_per_trial),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
    drgViolinPlot([1 2],[Cxy_pre_per_trial;Cxy_post_per_trial],edges,rand_offset,'k','k',4,1);
    
    title('Theta imaginary coherence')
    xticks([1 2])
    xticklabels({'Pre', 'Odor'})
    ylim([-0.5 1])
               
    %Enter the data for t-test/ranksum
    ii_rank=ii_rank+1;
    input_data(ii_rank).data=Cxy_post_per_trial;
    input_data(ii_rank).description='post odor';
    
    fprintf(1, ['\n\nRanksum or t-test for imginary coherence per trial\n'])
    %Now do the ranksums
    output_data = drgMutiRanksumorTtest(input_data);
    
end

pffft=1

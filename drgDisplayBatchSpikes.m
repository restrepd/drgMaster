function drgDisplayBatchSpikes(handles)

%This analyzes spike responses 
%Note: This not a general purpose program. Currently it analyzes Liz's data

close all
warning('off')


%Ask user for the drgb output .mat file and load those data
[handles.drgb.outFileName,handles.PathName] = uigetfile('*.mat','Select the drgb output file');
load([handles.PathName handles.drgb.outFileName])


%THESE VALUES ARE IMPORTANT
%VERY IMPORTANT: This is the index for this event in your handles.drgbchoices.evTypeNos
eventType=[2 5];
no_event_types=2;
%eventType=1;
%no_event_types=1;

evTypeLabels={'Hit';'CR'};
%evTypeLabels={'tstart'};

%Which percent correct bins do you want to use?
%For Liz's data we will only do 80-100, expert
percent_low=80;
percent_high=100;

%Pre time window
pre_start=-0.3;
pre_end=0.1;
pre_times=(handles_drgb.drgb.time>=pre_start)&(handles_drgb.drgb.time<=pre_end);

%Odor window
post_start=0.1;
post_end=0.4;
post_times=(handles_drgb.drgb.time>=post_start)&(handles_drgb.drgb.time<=post_end);

 
%First find ot how many odor-unit pairs responded
groupNo_oup=[];
no_oup=0;
for unitNo=1:handles_drgb.drgb.unit_no
    for evTy=1:no_event_types
        
        these_trials=(handles_drgb.drgb.unit(unitNo).perCorr>=80)&(handles_drgb.drgb.unit(unitNo).which_event(eventType(evTy),:)==1);
        if sum(these_trials)>3
            FR_pre=mean(handles_drgb.drgb.unit(unitNo).PSTH_per_trial(these_trials,pre_times),2);
            FR_post=mean(handles_drgb.drgb.unit(unitNo).PSTH_per_trial(these_trials,post_times),2);
            no_oup=no_oup+1;
            p_val_oup(no_oup)=ranksum(FR_pre,FR_post);
            groupNo_oup(no_oup)=handles_drgb.drgb.unit(unitNo).groupNo;
            normPost(no_oup)=mean(FR_post./FR_pre);
            perViol(no_oup)=handles_drgb.drgb.unit(unitNo).perViol;
        end
        
    end
end

pFDR=drsFDRpval(p_val_oup);

per_resp_WT=100*sum(p_val_oup(groupNo_oup==1)<=pFDR)/sum(groupNo_oup==1);
per_resp_Null=100*sum(p_val_oup(groupNo_oup==2)<=pFDR)/sum(groupNo_oup==2);
[p_WT_Null_oup, Q]= chi2test([sum(p_val_oup(groupNo_oup==1)<=pFDR), sum(groupNo_oup==1)-sum(p_val_oup(groupNo_oup==1)<=pFDR);...
    sum(p_val_oup(groupNo_oup==2)<=pFDR), sum(groupNo_oup==2)-sum(p_val_oup(groupNo_oup==2)<=pFDR)])

%Generate a bar graph
figure(1)
percent_WT=100*sum(p_val_oup(groupNo_oup==1)<=pFDR)/sum(groupNo_oup==1);
percent_Null=100*sum(p_val_oup(groupNo_oup==2)<=pFDR)/sum(groupNo_oup==2);
bar(1, percent_WT, 'b')
hold on
bar(2, percent_Null, 'r')
xlim([0.5 2.5])
title('Percent responding to odor, blue: WT, red: Null')
 
%Now do a pseudocolor plot of the responses
no_oup=0;
norm_oup=0;
for unitNo=1:handles_drgb.drgb.unit_no
    for evTy=1:no_event_types
        
        these_trials=(handles_drgb.drgb.unit(unitNo).perCorr>=80)&(handles_drgb.drgb.unit(unitNo).which_event(eventType(evTy),:)==1);
        if (sum(these_trials)>3)
            no_oup=no_oup+1;
            if p_val_oup(no_oup)<=pFDR
                norm_oup=norm_oup+1;
                normPSTH(norm_oup,:)=mean(handles_drgb.drgb.unit(unitNo).PSTH_per_trial(these_trials,:),1)/mean(mean(handles_drgb.drgb.unit(unitNo).PSTH_per_trial(these_trials,pre_times),2));
                groupNo_Norm_oup(norm_oup)=handles_drgb.drgb.unit(unitNo).groupNo;
                normResponse(norm_oup)=mean(normPSTH(norm_oup,post_times),2);
                normPost_sig(norm_oup)=normPost(no_oup);
            end
        end
        
    end
end

%Sort the unit odor pairs according to their pre
figure(2)
gr1_pre=normResponse(groupNo_Norm_oup==1)';
no_gr1=sum(groupNo_Norm_oup==1);
uops=[1:no_gr1]';
to_sort=[uops gr1_pre];
sorted_mtx=sortrows(to_sort,2);
PSTHgr1=normPSTH(groupNo_Norm_oup==1,:);
for ii=1:no_gr1
    jj=sorted_mtx(ii,1);
    sortedPSTHgr1(ii,:)=PSTHgr1(jj,:);
end
no_gr1=sum(groupNo_Norm_oup==1);
drg_pcolor(repmat(handles_drgb.drgb.time,no_gr1,1),repmat([1:no_gr1],length(handles_drgb.drgb.time),1)',sortedPSTHgr1)
colormap hot
caxis([0 3.3]);
title('Wild type reposnes for unit odor pairs')
xlabel('Time (sec)')

figure(3)
gr2_pre=normResponse(groupNo_Norm_oup==2)';
no_gr2=sum(groupNo_Norm_oup==2);
uops=[1:no_gr2]';
to_sort=[uops gr2_pre];
sorted_mtx=sortrows(to_sort,2);
PSTHgr2=normPSTH(groupNo_Norm_oup==2,:);
for ii=1:no_gr2
    jj=sorted_mtx(ii,1);
    sortedPSTHgr2(ii,:)=PSTHgr2(jj,:);
end
no_gr2=sum(groupNo_Norm_oup==2);
drg_pcolor(repmat(handles_drgb.drgb.time,no_gr2,1),repmat([1:no_gr2],length(handles_drgb.drgb.time),1)',sortedPSTHgr2)
colormap hot
caxis([0 3.3]);
title('Null reposnes for unit odor pairs')



pffft=1
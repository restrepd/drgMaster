function drgEventRelatedROCAnalysis(handles)

%Performs a ROC analysis of the power at time of the event
evTypeNo=handles.evTypeNo;
display_data=handles.displayData;
handles.displayData=0;

%Get the ERp for the first event
[log_P_t,no_trials_w_event,which_event,f,out_times,times,phase_per_trial,no_trials,no_events_per_trial,t_per_event_per_trial]=drgEventRelatedAnalysis(handles);
lpt1=zeros(no_trials_w_event,length(f));
lpt1(:,:)=log_P_t(no_events_per_trial>0,:,floor(length(out_times)/2)+1);
% lpt1(:,:)=mean(log_P_t(:,:,:),3);
%lpt1(:,:)=log_P_t(:,:,end);
meanlpt1=mean(lpt1(:,(f>=handles.burstLowF)&(f<=handles.burstHighF)),2);

%Get power for evTypeNo2
evTypeNo1=handles.evTypeNo;
handles.evTypeNo=handles.evTypeNo2;


%Get the ERp for the first event
[log_P_t,no_trials_w_event,which_event,f,out_times,times,phase_per_trial,no_trials,no_events_per_trial,t_per_event_per_trial]=drgEventRelatedAnalysis(handles);
lpt2=zeros(no_trials_w_event,length(f));
lpt2(:,:)=log_P_t(no_events_per_trial>0,:,floor(length(out_times)/2)+1);
% lpt2(:,:)=mean(log_P_t(:,:,:),3);
%lpt2(:,:)=log_P_t(:,:,end);
meanlpt2=mean(lpt2(:,(f>=handles.burstLowF)&(f<=handles.burstHighF)),2);

handles.evTypeNo=evTypeNo;
handles.displayData=display_data;

%Do the ROC analysis
roc_data=[meanlpt1; meanlpt2];
roc_data(1:length(roc_data),2)=[zeros(length(meanlpt1),1); ones(length(meanlpt2),1)];

%Calculate roc
roc_out=roc_calc(roc_data);

pffft=1;

end


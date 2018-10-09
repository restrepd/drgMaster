function drgLFProc(handles)



time_start=handles.time_start;
time_end=handles.time_end;
tstart=time_start+handles.time_pad;
tend=time_end-handles.time_pad;
tstart_ref=handles.startRef+handles.time_pad;
tend_ref=handles.endRef-handles.time_pad;
handles.time_start=handles.startRef-handles.window/2;
handles.time_end=handles.time_end+handles.window/2;
                

%Get power for evTypeNo
[t1,f1,all_Power1,all_Power_ref1, all_Power_timecourse1, this_trialNo1, perCorr1,which_event1]=drgGetLFPPowerForThisEvTypeNo(handles);
sz_all_power1=size(all_Power1);
this_all_power1=zeros(sz_all_power1(1),sz_all_power1(2));
this_all_power1(:,:)=mean(all_Power_timecourse1(:,:,(t1>=tstart)&(t1<=tend)),3);
this_all_dBpower1=10*log10(this_all_power1);
this_ref_power1=zeros(sz_all_power1(1),sz_all_power1(2));
this_ref_power1(:,:)=mean(all_Power_timecourse1(:,:,(t1>=tstart_ref)&(t1<=tend_ref)),3);
this_ref_dBpower1=10*log10(this_ref_power1);
 
%Get power for evTypeNo2
evTypeNo1=handles.evTypeNo;
evTypeNo2=handles.evTypeNo2;
handles.evTypeNo=handles.evTypeNo2;

[t2,f2,all_Power2,all_Power_ref2, all_Power_timecourse2, this_trialNo2, perCorr2,which_event2]=drgGetLFPPowerForThisEvTypeNo(handles);
sz_all_power2=size(all_Power2);
this_all_power2=zeros(sz_all_power2(1),sz_all_power2(2));
this_all_power2(:,:)=mean(all_Power_timecourse2(:,:,(t2>=tstart)&(t2<=tend)),3);
this_all_dBpower2=10*log10(this_all_power2);
this_ref_power2=zeros(sz_all_power2(1),sz_all_power2(2));
this_ref_power2(:,:)=mean(all_Power_timecourse2(:,:,(t2>=tstart_ref)&(t2<=tend_ref)),3);
this_ref_dBpower2=10*log10(this_ref_power2);
                
handles.evTypeNo=evTypeNo1;

roc_data=[];

size_allp1=size(all_Power1);
if handles.subtractRef==1
   roc_data(1:size_allp1(1),1)= mean(this_all_dBpower1-this_ref_dBpower1,2);
else
   roc_data(1:size_allp1(1),1)= mean(this_all_dBpower1,2);
end
roc_data(1:size_allp1(1),2)=zeros(size_allp1(1),1);

size_allp2=size(all_Power2);
if handles.subtractRef==1
   roc_data(size_allp1(1)+1:size_allp1(1)+size_allp2(1),1)= mean(this_all_dBpower2-this_ref_dBpower2,2);
else
   roc_data(size_allp1(1)+1:size_allp1(1)+size_allp2(1),1)= mean(this_all_dBpower2,2);
end
roc_data(size_allp1(1)+1:size_allp1(1)+size_allp2(1),2)=ones(size_allp2(1),1);

%Plot the histograms of delta dB
%Plot the lick frequency
try
    close 1
catch
end

hFig1 = figure(1);
set(hFig1, 'units','normalized','position',[.07 .1 .75 .3])
hold on

min_edge=min(roc_data(:,1))-0.1*(max(roc_data(:,1))-min(roc_data(:,1)));
max_edge=max(roc_data(:,1))+0.1*(max(roc_data(:,1))-min(roc_data(:,1)));
edges=[min_edge:(max_edge-min_edge)/20:max_edge];

h2=histogram(roc_data(size_allp1(1)+1:size_allp1(1)+size_allp2(1),1),edges);
h1=histogram(roc_data(1:size_allp1(1),1),edges);
legend(handles.drg.draq_d.eventlabels{evTypeNo2},handles.drg.draq_d.eventlabels{evTypeNo1})
ylim([0 1.2*max([max(h2.Values) max(h1.Values)])])
xlabel('Delta power (dB)')
ylabel('No of trials')
title('Histograms for delta power pre trial')


%Calculate roc
roc_out=roc_calc(roc_data);

handles.time_start=time_start;
handles.time_end=time_end;



fprintf(1, ['drgLFProc run from %d Hz to %d Hz\n\n'],handles.burstLowF,handles.burstHighF );



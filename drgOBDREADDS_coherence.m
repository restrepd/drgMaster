function handles_out=drgOBDREADDS_coherence(handles_choices)
%drgOBChr2LFPpower 
%provides power analysis for Joe's ChR2 runs
 
if exist('handles_choices')==0
    clear all
    close all

    %     handles.peakLFPNo=9;%This is the LFP number that will be processed
    %     electrode_label='Right OB';

%     handles.peakLFPNo=1;%This is the LFP number that will be processed
%     electrode_label='Right CA1';

    handles.peakLFPNo=9;%This is the LFP1 that will be processed
    handles.burstLFPNo=16;%This is the LFP2 that will be processed

    %1 Right hippocampus
    %8 Left hippocampus
    %9 Right olfactory bulb
    %16 Left olfactory bulb


    electrode_label='Left CA1';

    prefix_save_file='coh_hgt'; %Always use prp_

    % handles.peakLFPNo=[1 8 9 16];
    show_f_bandwidth=[65 95];

    %Choices for PAC frequency bandwidth
    handles.peakLowF=6;         %make peak and burst frequencies the same, usually 6-14 Hz (theta) 
    handles.peakHighF=14;

    handles.burstLowF=handles.peakLowF;       %burst is the power bandwidth (e.g. 65-95 Hz is high gamma)
    handles.burstHighF=handles.peakHighF;

    handles.n_phase_bins=50;

    handles.which_method=1;
    handles.save_drgb=0;
    handles.use_peakAngle=0;

    handles.window=1; %This is the FFT window in sec
    handles.noverlap=handles.window*0.9;

    %Load file
    % jtPathNames{1}='/Users/restrepd/Documents/Projects/Joe_OB_to_hippo/CNO/20231213_Siegfried_C5957_DREADD_n2_NewCNO_Pre_During_2mg_kg_2-Undecanone/20231213_Siegfried_C5957_DREADD_n2_NewCNO_Pre_231213_105533/';
    % jtFileNames{1}='jt_times_20231213_Siegfried_C5957_DREADD_n2_NewCNO_Pre_231213_105533.mat';
    % jtPathNames{2}='/Users/restrepd/Documents/Projects/Joe_OB_to_hippo/CNO/20231213_Siegfried_C5957_DREADD_n2_NewCNO_Pre_During_2mg_kg_2-Undecanone/20231213_Siegfried_C5957_DREADD_n2_NewCNO_2mg_kg_2-Undecanone_231213_111424/';
    % jtFileNames{2}='jt_times_20231213_Siegfried_C5957_DREADD_n2_NewCNO_2mg_kg_2-Undecanone_231213_111424.mat';
    %

    % jtPathNames{1}='/Volumes/Diego HD/Joe/Optogenetics/5xFADvsWT_1_hour_treatmentVsnone/Curley_WT_Treated/20240206_WT_Curley_Tx_2/20240206_WT_Curley_Tx_2_240206_120854/';
    % jtFileNames{1}='jt_times_20240206_WT_Curley_Tx_2_240206_120854.mat';


    jtPathNames{1}='/Volumes/Diego HD/Joe/Optogenetics/5xFADvsWT_1_hour_treatmentVsnone/Lennie_5xFAD_Treated/20240206_5xFAD_Lennie_Tx_2/20240206_5xFAD_Lennie_Tx_2_240206_133005/';
    jtFileNames{1}='jt_times_20240206_5xFAD_Lennie_Tx_2_240206_133005.mat';

    %
    %     [jtFileName,jtPathName] = uigetfile('jt_times*.mat','Select jt_times file to open');
%     handles.jtfullName=[jtPathName,jtFileName];
    handles.jtFileNames=jtFileNames;
    handles.jtPathNames=jtPathNames;

    handles.showData=1;
    

    ii_laser_start=60*5;
    ii_laser_end=60*65;

else
    %Use this if the file is called by another function
    show_f_bandwidth=[35 45];
    jtFileName=handles_choices.jtFileName;
    jtPathName=handles_choices.jtPathName;
    electrode_label=handles_choices.electrode_label;
    handles.peakLFPNo=handles_choices.peakLFPNo;%This is the LFP number that will be processed
    % handles.peakLFPNo=[1 8 9 16];
 
   
    handles.burstLowF=handles_choices.burstLowF;
    handles.burstHighF=handles_choices.burstHighF;
    handles.showData=0;

        %Load file

%     handles.jtfullName=[jtPathName,jtFileName];
    handles.jtFileNames=jtFileNames;
    handles.jtPathNames=jtPathNames;
  
end

handles_out=[];

handles.displayData=1;


% try
%     mkdir('figures')
% catch
% end

%Initialize handles
handles.sessionNo=1;
sessionNo=1;
handles.unitNo=1;
handles.evTypeNo=1; %This is event 1 in continuous
handles.notch60=1;
handles.subtractRef=0; %Subtract the mean reference power?
handles.data_vs_simulate=0;

handles.autoscale=1;





% handles.burstLFPNo=handles.peakLFPNo;

%Note: for contnuous
% t_start is 3.2 sec after the start of the trial and
% trials are spaced by 9 sec
%These times are set so that the entire trial is processed
%These settings read 9 sec
handles.time_start=-3.2-0.2;
handles.time_end=(8-3.2)+0.2;

handles.delta_trial=9;

%This reference is not used here
handles.startRef=-2.2;
handles.endRef=0.2;

handles.time_pad=0.2;

handles.dt_tPRP=0.03333;


%Initialize variables
digital_LFPNo=22;
dec_camera=[];
dec_laser=[];
bit_camera=[];
bit_laser=[];
max_time_bins=120*60;

%Decimate the power further to 1 per sec
mean_window=1; %window to average in seconds

ii_t=0;

handles_out=[];
handles_out.coherence=[];
handles_out.img_coherence=[];

time=[];
tNum_all=0;
ii_start_per_file=[];
ii_end_per_file=[];
for fileNo=1:length(jtFileNames)

    jtPathName=jtPathNames{fileNo};
    jtFileName=jtFileNames{fileNo};
    cd(jtPathName(1:end-1))

    drgRead_jt_times(jtPathName,jtFileName)
    FileName=[jtFileName(10:end-4) '_drg.mat'];
    handles.fullName=[jtPathName,FileName];
    load(handles.fullName);
    handles.drg=drg;
    handles.drg.drta_p.fullName=[jtPathName,handles.drg.drta_p.FileName];

    handles.trialNo=1;
    handles.lastTrialNo=handles.drg.drta_p.trialNo;
    dec_n=fix(handles.drg.session(sessionNo).draq_p.ActualRate/1000);

    % handles=drgLFPwaveSpectrogramRefOtherTrial(handles);

    %First figure out when the laser input and camera coverage takes place

    % for trNo=1:handles.lastTrialNo
    %     [digital_input, trialNo, can_read] = drgGetTrialLFPData(handles, digital_LFPNo, trNo, handles.evTypeNo, handles.time_start, handles.time_end);
    %     bit_camera=bitget(uint16(digital_input), 1, 'uint16');
    %     bit_laser=bitget(uint16(digital_input), 2, 'uint16');
    %     dec_camera=[dec_camera decimate(double(bit_camera(1:180000)),dec_n)];
    %     dec_laser=[dec_laser decimate(double(bit_laser(1:180000)),dec_n)];
    % end
    all_Cxy_timecourse=[];
    [t,f, all_Cxy_timecourse, this_trialNo]=drgGetLFPCoherenceForThisEvTypeNo(handles);

    this_per_trialCxy=zeros(1,size(all_Cxy_timecourse,1));
    this_per_trialCxy(1,:)=mean(mean(all_Cxy_timecourse,3),2)';

    handles_out.coherence=[handles_out.coherence this_per_trialCxy];

    img_all_Cxy_timecourse=[];
    [t,f, img_all_Cxy_timecourse, this_trialNo]=drgGetLFPimgCoherenceForThisEvTypeNo(handles);

    this_per_trialCxyimg=zeros(1,size(img_all_Cxy_timecourse,1));
    this_per_trialCxyimg(1,:)=mean(mean(img_all_Cxy_timecourse,3),2)';

    handles_out.img_coherence=[handles_out.img_coherence this_per_trialCxyimg];

    if fileNo==1
        ii_start_per_file(fileNo)=1;
    else
        ii_start_per_file(fileNo)=ii_end_per_file(fileNo-1)+1;
    end
    ii_end_per_file(fileNo)=length(handles_out.img_coherence);

    for trNum=1:length(this_trialNo)

        tNum_all=tNum_all+1;


        if tNum_all==1
            time(1)=0;
        else
            time(tNum_all)=time(tNum_all-1)+handles.delta_trial;
        end

        if (trNum==1)&(fileNo==2)
            time(tNum_all)=time(tNum_all)+handles.delta_trial*4*60/9;
        end
    end

    % ii_end_per_file(fileNo)=tNum_all;
end

% if fileNo==1
%     log_P_timecourse=zeros(length(f),max_time_bins);
% end

% for trNo=1:handles.lastTrialNo
%     for ii_window=1:size(all_Power_timecourse,3)/1000
%         ii_t=ii_t+1;
%         this_mean_logP=zeros(length(f),1);
%         ii_from=(ii_window-1)*mean_window*1000+1;
%         ii_to=ii_window*mean_window*1000;
%         this_mean_logP(:,1)=mean(10*log10(all_Power_timecourse(trNo,:,ii_from:ii_to)),3);
%         log_P_timecourse(:,ii_t)=this_mean_logP;
%     end
% end
% if fileNo==1
%     ii_laser_start=ii_t;
% end
% ii_t=ii_t+180;
%  if fileNo==1
%     ii_laser_end=ii_t;
% end

%
% ii_laser_start=floor(find(dec_laser>0.5,1,'first')/(1000*mean_window));
% ii_laser_end=ceil(find(dec_laser>0.5,1,'last')/(1000*mean_window));

handles_out.time=time;
handles_out.ii_start_per_file=ii_start_per_file;
handles_out.ii_end_per_file=ii_end_per_file;

save_base=jtFileNames{1};
save([jtPathNames{1} prefix_save_file '_ch' num2str(handles.peakLFPNo), '_', save_base(9:end)],'handles','handles_out')

close all
figNo=0;

%Plot coherence
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig2 = figure(figNo);
set(hFig2, 'units','normalized','position',[.25 .1 .4 .25])
hold on

coherence=handles_out.coherence;

for fileNo=1:length(jtFileNames)
    plot(time(ii_start_per_file(fileNo):ii_end_per_file(fileNo)),coherence(ii_start_per_file(fileNo):ii_end_per_file(fileNo)),'ob-')
end


if length(coherence)>20
    no_conv_points=6;
    conv_win=ones(1,no_conv_points);
    coherence_extend=[mean(coherence(1:no_conv_points/2))*ones(1,no_conv_points) coherence mean(coherence(end-(no_conv_points/2)+1:end))*ones(1,no_conv_points)];
    conv_coherence_extend=conv(coherence_extend,conv_win)/no_conv_points;
    conv_coherence=conv_coherence_extend(no_conv_points+1:no_conv_points+length(coherence));
    for fileNo=1:length(jtFileNames)
        plot(time(ii_start_per_file(fileNo):ii_end_per_file(fileNo)),conv_coherence(ii_start_per_file(fileNo):ii_end_per_file(fileNo)),'-b','LineWidth',3)
    end

end
theseyl=ylim;
% between_ii=(time(ii_end_per_file(1))+time(ii_start_per_file(length(jtFileNames))))/2;
% plot([between_ii between_ii],theseyl,'-k','LineWidth',3)
xlabel('Time (sec)')
ylabel('Coherence')
title('Coherence vs time')

%Plot imaginary coherence
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig2 = figure(figNo);
set(hFig2, 'units','normalized','position',[.25 .1 .4 .25])
hold on

img_coherence=handles_out.img_coherence;

for fileNo=1:length(jtFileNames)
    plot(time(ii_start_per_file(fileNo):ii_end_per_file(fileNo)),img_coherence(ii_start_per_file(fileNo):ii_end_per_file(fileNo)),'ob-')
end


if length(img_coherence)>20
    no_conv_points=6;
    conv_win=ones(1,no_conv_points);
    img_coherence_extend=[mean(img_coherence(1:no_conv_points/2))*ones(1,no_conv_points) img_coherence mean(img_coherence(end-(no_conv_points/2)+1:end))*ones(1,no_conv_points)];
    conv_img_coherence_extend=conv(img_coherence_extend,conv_win)/no_conv_points;
    conv_img_coherence=conv_img_coherence_extend(no_conv_points+1:no_conv_points+length(img_coherence));
    for fileNo=1:length(jtFileNames)
        plot(time(ii_start_per_file(fileNo):ii_end_per_file(fileNo)),conv_img_coherence(ii_start_per_file(fileNo):ii_end_per_file(fileNo)),'-b','LineWidth',3)
    end
end
theseyl=ylim;
% between_ii=(time(ii_end_per_file(1))+time(ii_start_per_file(length(jtFileNames))))/2;
% plot([between_ii between_ii],theseyl,'-k','LineWidth',3)
xlabel('Time (sec)')
ylabel('Imaginary Coherence')
title('Imaginary Coherence vs time')


pffft=1
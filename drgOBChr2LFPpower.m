
%drgOBChr2LFPpower 
close all
clear all

handles.peakLFPNo=8;%This is the LFP number that will be processed
show_f_bandwidth=[35 45];


electrode_labels{1}='Right CA1';
electrode_labels{8}='Left CA1';
electrode_labels{9}='Right OB';
electrode_labels{16}='Left OB';


%Initialize handles
handles.sessionNo=1;
sessionNo=1;
handles.unitNo=1;
handles.evTypeNo=1; %This is event 1 in continuous
handles.notch60=1;
handles.subtractRef=1; %Subtract the mean reference power?
handles.data_vs_simulate=0;
handles.displayData=1;
handles.autoscale=1;
 
handles.burstLowF=1;
handles.burstHighF=100;



handles.burstLFPNo=handles.peakLFPNo;

%Note: for contnuous
% t_start is 3.2 sec after the start of the trial and
% trials are spaced by 9 sec
%These times are set so that the entire trial is processed
%These settings read 9 sec
handles.time_start=-3.2-0.2;
handles.time_end=(9-3.2)+0.2;

%This reference is not used here
handles.startRef=-2.2;
handles.endRef=0.2;

handles.time_pad=0.2;





%Load file
[jtFileName,jtPathName] = uigetfile('jt_times*.mat','Select jt_times file to open');
handles.jtfullName=[jtPathName,jtFileName];
handles.jtFileName=jtFileName;
handles.jtPathName=jtPathName;
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
digital_LFPNo=22;
dec_camera=[];
dec_laser=[];
bit_camera=[];
bit_laser=[];
for trNo=1:handles.lastTrialNo
    [digital_input, trialNo, can_read] = drgGetTrialLFPData(handles, digital_LFPNo, trNo, handles.evTypeNo, handles.time_start, handles.time_end);
    bit_camera=bitget(uint16(digital_input), 1, 'uint16');
    bit_laser=bitget(uint16(digital_input), 2, 'uint16');
    dec_camera=[dec_camera decimate(double(bit_camera(1:180000)),dec_n)];
    dec_laser=[dec_laser decimate(double(bit_laser(1:180000)),dec_n)];
end



[out_t,f,all_Power, all_Power_ref, all_Power_timecourse, this_trialNo, perCorr_pertr, which_event, is_not_moving]=drgGetLFPwavePowerForThisEvTypeNo(handles);

%Now decimate the power further to 1 per sec
mean_window=1; %window to average in seconds

log_P_timecourse=zeros(length(f),handles.lastTrialNo*size(all_Power_timecourse,3)/(mean_window*1000));
ii_t=0;
for trNo=1:handles.lastTrialNo
    for ii_window=1:size(all_Power_timecourse,3)/1000
        ii_t=ii_t+1;
        this_mean_logP=zeros(length(f),1);
        ii_from=(ii_window-1)*mean_window*1000+1;
        ii_to=ii_window*mean_window*1000;
        this_mean_logP(:,1)=mean(10*log10(all_Power_timecourse(trNo,:,ii_from:ii_to)),3);
        log_P_timecourse(:,ii_t)=this_mean_logP;
    end
end

ii_laser_start=floor(find(dec_laser>0.5,1,'first')/(1000*mean_window));
ii_laser_end=ceil(find(dec_laser>0.5,1,'last')/(1000*mean_window));

%Subtract reference
this_mean_reference_logP=zeros(length(f),1);
this_mean_reference_logP(:,1)=mean(log_P_timecourse(:,1:ii_laser_start-1),2);
log_P_timecourse_ref=repmat(this_mean_reference_logP,1,size(log_P_timecourse,2));
log_P_timecourse=log_P_timecourse-log_P_timecourse_ref;
time=[mean_window:mean_window:mean_window*size(log_P_timecourse,2)];


if handles.displayData==1
    
    if handles.autoscale==1
        maxLogPper=prctile(log_P_timecourse(:),99);
        minLogPper=prctile(log_P_timecourse(:),1);
        %Note: Diego added this on purpose to limit the range to 10 dB
        %This results in emphasizing changes in the top 10 dB
        if maxLogPper-minLogPper>16
            minLogPper=maxLogPper-16;
        end
    else
        maxLogPper=handles.maxLogP;
        minLogPper=handles.minLogP;
    end

    figNo=0;

    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    %Plot the timecourse
    hFig = figure(figNo);
    set(hFig, 'units','normalized','position',[.07 .5 .75 .3])

    drg_pcolor(repmat(time,length(f),1)',repmat(f,length(time),1),log_P_timecourse')

    colormap fire
    shading interp
    caxis([minLogPper maxLogPper]);

    hold on 
    plot([time(ii_laser_start) time(ii_laser_start)],[f(1) f(end)],'-k','LineWidth',2)
    plot([time(ii_laser_end) time(ii_laser_end)],[f(1) f(end)],'-k','LineWidth',2)

    xlabel('Time (sec)')
    ylabel('Frequency (Hz)');
    title(['Wavelet power (dB) ' electrode_labels{handles.peakLFPNo} ' ' jtFileName])


    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);
    set(hFig, 'units','normalized','position',[.83 .5 .05 .3])

    prain=[minLogPper:(maxLogPper-minLogPper)/99:maxLogPper];
    drg_pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
    colormap fire
    shading interp
    ax=gca;
    set(ax,'XTickLabel','')

    %Plot the dependence on frequency before, during and after laser
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    hFig = figure(figNo);
    set(hFig, 'units','normalized','position',[.07 .15 .25 .25])
    hold on

    CI=[];
    CI = bootci(1000, {@mean, log_P_timecourse(:,1:ii_laser_start)'})';
    CI(:,1)= mean(log_P_timecourse(:,1:ii_laser_start)')'-CI(:,1);
    CI(:,2)=CI(:,2)- mean(log_P_timecourse(:,1:ii_laser_start)')';
    [hlCR, hpCR] = boundedline(f',mean(log_P_timecourse(:,1:ii_laser_start)'), CI, 'b');

    CI=[];
    CI = bootci(1000, {@mean, log_P_timecourse(:,ii_laser_end:end)'})';
    CI(:,1)= mean(log_P_timecourse(:,ii_laser_end:end)')'-CI(:,1);
    CI(:,2)=CI(:,2)- mean(log_P_timecourse(:,ii_laser_end:end)')';
    [hlCR, hpCR] = boundedline(f',mean(log_P_timecourse(:,ii_laser_end:end)')', CI, 'c');

    CI=[];
    CI = bootci(1000, {@mean, log_P_timecourse(:,ii_laser_start:ii_laser_end)'})';
    CI(:,1)= mean(log_P_timecourse(:,ii_laser_start:ii_laser_end)')'-CI(:,1);
    CI(:,2)=CI(:,2)- mean(log_P_timecourse(:,ii_laser_start:ii_laser_end)')';
    [hlCR, hpCR] = boundedline(f',mean(log_P_timecourse(:,ii_laser_start:ii_laser_end)')', CI, 'm');

    title('Wavelet power spectrum, blue=before, magenta=laser on, cyan=after')
    xlabel('Frequency (Hz)')
    ylabel('dB')

    %Now show the timecourse for the show frequency bandwidth
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    hFig = figure(figNo);
    set(hFig, 'units','normalized','position',[.37 .15 .25 .25])
    hold on

    this_log_P_timecourse=zeros(1,size(log_P_timecourse,2));
    this_log_P_timecourse(1,:)=mean(log_P_timecourse((f>=show_f_bandwidth(1))&(f<=show_f_bandwidth(2)),:));
    plot(time, this_log_P_timecourse,'-k')

    this_ylim=ylim;
    hold on
    plot([time(ii_laser_start) time(ii_laser_start)],[this_ylim(1) this_ylim(2)],'-k','LineWidth',2)
    plot([time(ii_laser_end) time(ii_laser_end)],[this_ylim(1) this_ylim(2)],'-k','LineWidth',2)

    title(['Wavelet power timecourse for bandwidth from ' num2str(show_f_bandwidth(1)) ' to ' num2str(show_f_bandwidth(2)) ' Hz'])
    xlabel('Time (sec)')
    ylabel('dB')
   

end


function handles_out=drgOBtoHIP_LFPpowerv2(handles)
%drgOBChr2LFPpower 
%provides power analysis for Joe's ChR2 runs

if exist('handles')==0
    clear all
    close all

    %     handles.peakLFPNo=9;%This is the LFP number that will be processed
    %     electrode_label='Right OB';

    %     handles.peakLFPNo=1;%This is the LFP number that will be processed
    %     electrode_label='Right CA1';

    handles.peakLFPNo=9;%This is the LFP number that will be processed
    handles.electrode_label{1}='Right HIP';
    handles.electrode_label{8}='Left HIP';
    handles.electrode_label{9}='Right OB';
    handles.electrode_label{16}='Left OB';
    electrode_label=handles.electrode_label{handles.peakLFPNo};

    %1 Right hippocampus
    %8 Left hippocampus
    %9 Right olfactory bulb
    %16 Left olfactory bulb

    % handles.peakLFPNo=[1 8 9 16];
   


    handles.window=1; %This is the FFT window in sec  %Tort is 1 sec, old DR 0.37
    handles.noverlap=handles.window*0.9; 
  
    handles.burstLowF=1;
    handles.burstHighF=100;

    %Load file
    % handles.jtPathNames{1}='/Users/restrepd/Documents/Projects/Joe_OB_to_hippo/CNO/20231213_Siegfried_C5957_DREADD_n2_NewCNO_Pre_During_2mg_kg_2-Undecanone/20231213_Siegfried_C5957_DREADD_n2_NewCNO_Pre_231213_105533/'
    % handles.jtFileNames{1}='jt_times_20231213_Siegfried_C5957_DREADD_n2_NewCNO_Pre_231213_105533.mat';
    % handles.jtPathNames{2}='/Users/restrepd/Documents/Projects/Joe_OB_to_hippo/CNO/20231213_Siegfried_C5957_DREADD_n2_NewCNO_Pre_During_2mg_kg_2-Undecanone/20231213_Siegfried_C5957_DREADD_n2_NewCNO_2mg_kg_2-Undecanone_231213_111424/'
    % handles.jtFileNames{2}='jt_times_20231213_Siegfried_C5957_DREADD_n2_NewCNO_2mg_kg_2-Undecanone_231213_111424.mat';
    % 
    % handles.ii_laser_start=770;
    % handles.ii_laser_end=5700;

    %Used for the grant ChR2 WT
    % handles.jtPathNames{1}='/Volumes/Diego HD/Joe/Optogenetics/5xFADvsWT_1_hour_treatmentVsnone/Curley_WT_Treated/20240206_WT_Curley_Tx_2/20240206_WT_Curley_Tx_2_240206_120854/';
    % handles.jtFileNames{1}='jt_times_20240206_WT_Curley_Tx_2_240206_120854.mat';

    % %Used for the grant no ChR2 WT
    % handles.jtPathNames{1}='/Volumes/Diego HD/Joe/022824_ChR2WT_withLaser_240228_085402/';
    % handles.jtFileNames{1}='jt_times_022824_ChR2WT_withLaser_240228_085402.mat';

    % handles.ii_laser_start=250;
    % handles.ii_laser_end=3500;

    %Joe closed loop
    % handles.jtPathNames{1}='/Users/restrepd/Documents/Projects/Closed loop/Joe_ClosedLoop/20240306_Donatello_Closed_Loop_Test_withBandPass_4_240306_112305/';
    % handles.jtFileNames{1}='jt_times_20240306_Donatello_Closed_Loop_Test_withBandPass_4_240306_112305.mat';

    %Lyra day 1
    % handles.jtPathNames{1}='/Volumes/Diego Mac Drive/5xFADvsWT_1_hour_treatmentVsnone/N_4/Lyra 5xFAD WT Tx/20240819_Lyra_Day1_Tx/20240819_Lyra_Day1_Tx_1_240819_094218/';
    % handles.jtFileNames{1}='jt_times_20240819_Lyra_Day1_Tx_1_240819_094218.mat';

    %Post week 2 202501117
    %Lyra day 1
    handles.jtPathNames{1}='/Users/restrepd/Documents/Projects/JoeV/Exosomes/20250117_Mouse4 VZV_250117_094111/';
    handles.jtFileNames{1}='jt_times_20250117_Mouse4 VZV_250117_094111.mat';

    %
    % handles.laser_start_time=5; %Start time in minutes
    % handles.laser_end_time=65; %Start time in minutes
    % 
    % handles.ii_laser_start=60*handles.laser_start_time;
    % handles.ii_laser_end=60*handles.laser_end_time;

    %If you want the program to ask you for files use this code
    % [jtFileName,jtPathName] = uigetfile('jt_times*.mat','Select jt_times file to open');
    % handles.jtfullName=[jtPathName,jtFileName];
    % 
    % handles.jtPathNames{1}=jtPathName;
    % handles.jtFileNames{1}=jtFileName;

    handles.showData=1;
  
end

handles_out=[];

show_f_bandwidth=[35 45];
before_f=[22 27];
after_f=[50 57];
dt_per_trial=9;
mean_window=1; %window to average in seconds

handles_out.show_f_bandwidth=show_f_bandwidth;
handles_out.before_f=before_f;
handles_out.after_f=after_f;
handles_out.dt_per_trial=dt_per_trial;
handles_out.mean_window=mean_window;

handles.displayData=1;


cd(handles.jtPathNames{1}(1:end-1))
try
    mkdir('figures')
catch
end

%Initialize handles
handles.sessionNo=1;
sessionNo=1;
handles.unitNo=1;
handles.evTypeNo=1; %This is event 1 in continuous
handles.notch60=1;
handles.subtractRef=0; %Subtract the mean reference power?
handles.data_vs_simulate=0;

handles.autoscale=1;
 




handles.burstLFPNo=handles.peakLFPNo;

%Note: for contnuous
% t_start is 3.2 sec after the start of the trial and
% trials are spaced by 9 sec
%These times are set so that the entire trial is processed
%These settings read 9 sec
handles.time_start=-3.2-0.2;
handles.time_end=(dt_per_trial-3.2)+0.2;

%This reference is not used here
handles.startRef=-2.2;
handles.endRef=0.2;

handles.time_pad=0.2;


%Initialize variables
digital_LFPNo=22;
dec_camera=[];
dec_laser=[];
bit_camera=[];
bit_laser=[];
max_time_bins=120*60;

%Decimate the power further to 1 per sec


ii_t=0;
file_starts=[];

for fileNo=1:length(handles.jtFileNames)

    jtPathName=handles.jtPathNames{fileNo};
    jtFileName=handles.jtFileNames{fileNo};
    jtfullName=[jtPathName,jtFileName];
    file_info = dir(jtfullName);
        if fileNo==1
        start_date=datenum(file_info.date);
    end
    file_starts(fileNo)=(datenum(file_info.date)-start_date)*24*60*60;
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

    for trNo=1:handles.lastTrialNo
        [digital_input, trialNo, can_read] = drgGetTrialLFPData(handles, digital_LFPNo, trNo, handles.evTypeNo, handles.time_start, handles.time_end);
        bit_camera=bitget(uint16(digital_input), 1, 'uint16');
        bit_laser=bitget(uint16(digital_input), 2, 'uint16');
        dec_camera=[dec_camera decimate(double(bit_camera(1:180000)),dec_n)];
        dec_laser=[dec_laser decimate(double(bit_laser(1:180000)),dec_n)];
    end



    [out_t,f,all_Power, all_Power_ref, all_Power_timecourse, this_trialNo, perCorr_pertr, which_event]=drgGetLFPwavePowerForThisEvTypeNo(handles);

    
    if fileNo==1
        log_P_timecourse=zeros(length(f),max_time_bins);
    end
    ii_t_file=file_starts(fileNo);
    ii_t=0;
    for trNo=1:handles.lastTrialNo
        for ii_window=1:dt_per_trial
            ii_t=ii_t+1;
            this_mean_logP=zeros(length(f),1);
            ii_from=(ii_window-1)*mean_window*1000+1;
            ii_to=ii_window*mean_window*1000;
            if ii_to>size(all_Power_timecourse,3)
                ii_to=size(all_Power_timecourse,3);
            end
            this_mean_logP(:,1)=mean(10*log10(all_Power_timecourse(trNo,:,ii_from:ii_to)),3);
            log_P_timecourse(:,ii_t+round(ii_t_file))=this_mean_logP;
        end
    end
    % handles.ii_laser_start=floor(find(dec_laser>0.5,1,'first')/(1000*mean_window));
    % if isempty(handles.ii_laser_start)
    %     handles.ii_laser_start=(min_laser_start/min_exp_end)*ceil(handles.lastTrialNo*size(all_Power_timecourse,3)/(mean_window*1000));
    % end
    % handles.ii_laser_end=ceil(find(dec_laser>0.5,1,'last')/(1000*mean_window));
    % if isempty(handles.ii_laser_end)
    %     handles.ii_laser_end=(min_laser_end/min_exp_end)*ceil(handles.lastTrialNo*size(all_Power_timecourse,3)/(mean_window*1000));
    % end
end
log_P_timecourse=log_P_timecourse(:,1:ii_t+round(ii_t_file));

%Subtract reference
% this_mean_reference_logP=zeros(length(f),1);
% this_mean_reference_logP(:,1)=mean(log_P_timecourse(:,1:handles.ii_laser_start-1),2);
% log_P_timecourse_ref=repmat(this_mean_reference_logP,1,size(log_P_timecourse,2));

% log_P_timecourse_no_sub=log_P_timecourse;
% log_P_timecourse=log_P_timecourse-log_P_timecourse_ref;
time=[mean_window:mean_window:mean_window*size(log_P_timecourse,2)]/60;

handles_out.log_P_timecourse=log_P_timecourse;
handles_out.time=time;
% handles_out.handles.ii_laser_start=handles.ii_laser_start;
% handles_out.handles.ii_laser_end=handles.ii_laser_end;
handles_out.f=f;

if handles.showData==1
    
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
%     plot([time(handles.ii_laser_start) time(handles.ii_laser_start)],[f(1) f(end)],'-k','LineWidth',2)
%     plot([time(handles.ii_laser_end) time(handles.ii_laser_end)],[f(1) f(end)],'-k','LineWidth',2)

    xlabel('Time (min)')
    ylabel('Frequency (Hz)');
    title(['Wavelet power (dB) ' handles.electrode_label{handles.peakLFPNo} ' ' jtFileName])

    fig_file_name=[jtPathName 'figures/spectrogram' num2str(handles.peakLFPNo) '.fig'];
    savefig(fig_file_name)
    
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

    % %Plot the dependence on frequency before, during and after laser
    % figNo=figNo+1;
    % try
    %     close(figNo)
    % catch
    % end
    % hFig = figure(figNo);
    % set(hFig, 'units','normalized','position',[.07 .15 .25 .25])
    % hold on
    % 
    % CI=[];
    % CI = bootci(1000, {@mean, log_P_timecourse(:,1:handles.ii_laser_start)'})';
    % CI(:,1)= mean(log_P_timecourse(:,1:handles.ii_laser_start)')'-CI(:,1);
    % CI(:,2)=CI(:,2)- mean(log_P_timecourse(:,1:handles.ii_laser_start)')';
    % [hlCR, hpCR] = boundedline(f',mean(log_P_timecourse(:,1:handles.ii_laser_start)'), CI, 'b');
    % 
    % % CI=[];
    % % CI = bootci(1000, {@mean, log_P_timecourse(:,handles.ii_laser_end:end)'})';
    % % CI(:,1)= mean(log_P_timecourse(:,handles.ii_laser_end:end)')'-CI(:,1);
    % % CI(:,2)=CI(:,2)- mean(log_P_timecourse(:,handles.ii_laser_end:end)')';
    % % [hlCR, hpCR] = boundedline(f',mean(log_P_timecourse(:,handles.ii_laser_end:end)')', CI, 'c');
    % 
    % CI=[];
    % CI = bootci(1000, {@mean, log_P_timecourse(:,handles.ii_laser_start:handles.ii_laser_end)'})';
    % CI(:,1)= mean(log_P_timecourse(:,handles.ii_laser_start:handles.ii_laser_end)')'-CI(:,1);
    % CI(:,2)=CI(:,2)- mean(log_P_timecourse(:,handles.ii_laser_start:handles.ii_laser_end)')';
    % [hlCR, hpCR] = boundedline(f',mean(log_P_timecourse(:,handles.ii_laser_start:handles.ii_laser_end)')', CI, 'm');
    % 
    % ylim([-3 5])
    % this_ylim=ylim;
    % 
    % title('Wavelet power spectrum, blue=before, magenta=laser on, cyan=after')
    % xlabel('Frequency (Hz)')
    % ylabel('dB')
    % 
    % %Plot the dependence on frequency before, during and after laser
    % figNo=figNo+1;
    % try
    %     close(figNo)
    % catch
    % end
    % hFig = figure(figNo);
    % set(hFig, 'units','normalized','position',[.07 .15 .25 .25])
    % hold on
    % 
    % CI=[];
    % CI = bootci(1000, {@mean, log_P_timecourse_no_sub(:,1:handles.ii_laser_start)'})';
    % CI(:,1)= mean(log_P_timecourse_no_sub(:,1:handles.ii_laser_start)')'-CI(:,1);
    % CI(:,2)=CI(:,2)- mean(log_P_timecourse_no_sub(:,1:handles.ii_laser_start)')';
    % [hlCR, hpCR] = boundedline(f',mean(log_P_timecourse_no_sub(:,1:handles.ii_laser_start)'), CI, 'b');
    % 
    % CI=[];
    % CI = bootci(1000, {@mean, log_P_timecourse_no_sub(:,handles.ii_laser_end:end)'})';
    % CI(:,1)= mean(log_P_timecourse_no_sub(:,handles.ii_laser_end:end)')'-CI(:,1);
    % CI(:,2)=CI(:,2)- mean(log_P_timecourse_no_sub(:,handles.ii_laser_end:end)')';
    % [hlCR, hpCR] = boundedline(f',mean(log_P_timecourse_no_sub(:,handles.ii_laser_end:end)')', CI, 'c');
    % 
    % CI=[];
    % CI = bootci(1000, {@mean, log_P_timecourse_no_sub(:,handles.ii_laser_start:handles.ii_laser_end)'})';
    % CI(:,1)= mean(log_P_timecourse_no_sub(:,handles.ii_laser_start:handles.ii_laser_end)')'-CI(:,1);
    % CI(:,2)=CI(:,2)- mean(log_P_timecourse_no_sub(:,handles.ii_laser_start:handles.ii_laser_end)')';
    % [hlCR, hpCR] = boundedline(f',mean(log_P_timecourse_no_sub(:,handles.ii_laser_start:handles.ii_laser_end)')', CI, 'm');
    % 
    % % ylim([-10 20])
    % % this_ylim=ylim;
    % 
    % title('Wavelet power spectrum (no subtraction), blue=before, magenta=laser on, cyan=after')
    % xlabel('Frequency (Hz)')
    % ylabel('dB')

    % %Now show the timecourse for the show frequency bandwidth
    % figNo=figNo+1;
    % try
    %     close(figNo)
    % catch
    % end
    % hFig = figure(figNo);
    % set(hFig, 'units','normalized','position',[.37 .15 .25 .25])
    % hold on
    % 
    % this_log_P_timecourse=zeros(1,size(log_P_timecourse,2));
    % this_log_P_timecourse(1,:)=mean(log_P_timecourse((f>=show_f_bandwidth(1))&(f<=show_f_bandwidth(2)),:));
    % plot(time, this_log_P_timecourse,'-k')
    % 
    % ylim([-10 20])
    % this_ylim=ylim;
    % hold on
    % plot([time(handles.ii_laser_start) time(handles.ii_laser_start)],[this_ylim(1) this_ylim(2)],'-k','LineWidth',2)
    % plot([time(handles.ii_laser_end) time(handles.ii_laser_end)],[this_ylim(1) this_ylim(2)],'-k','LineWidth',2)
    % 
    % title([handles.electrode_label{handles.peakLFPNo} ' Wavelet power timecourse for bandwidth from ' num2str(show_f_bandwidth(1)) ' to ' num2str(show_f_bandwidth(2)) ' Hz'])
    % xlabel('Time (min)')
    % ylabel('dB')
    % 
    % %Now show the timecourse for the 40 Hz peak
    % figNo=figNo+1;
    % try
    %     close(figNo)
    % catch
    % end
    % hFig = figure(figNo);
    % set(hFig, 'units','normalized','position',[.37 .15 .25 .25])
    % hold on
    % 
    % this_log_P_timecourse=zeros(1,size(log_P_timecourse,2));
    % this_log_P_timecourse(1,:)=mean(log_P_timecourse((f>=show_f_bandwidth(1))&(f<=show_f_bandwidth(2)),:));
    % f_b=mean(f((f>=show_f_bandwidth(1))&(f<=show_f_bandwidth(2))));
    % this_log_P_bef_timecourse=zeros(1,size(log_P_timecourse,2));
    % this_log_P_bef_timecourse(1,:)=mean(log_P_timecourse((f>=before_f(1))&(f<=before_f(2)),:));
    % bef_f=mean(f((f>=before_f(1))&(f<=before_f(2))));
    % this_log_P_af_timecourse=zeros(1,size(log_P_timecourse,2));
    % this_log_P_af_timecourse(1,:)=mean(log_P_timecourse((f>=after_f(1))&(f<=after_f(2)),:));
    % af_f=mean(f((f>=after_f(1))&(f<=after_f(2))));
    % 
    % 
    % slope=(this_log_P_af_timecourse-this_log_P_bef_timecourse)/(af_f-bef_f);
    % this_peak_log_P_af_timecourse=zeros(1,size(log_P_timecourse,2));
    % 
    % this_peak_log_P_af_timecourse=this_log_P_timecourse-(this_log_P_bef_timecourse+slope*(f_b-bef_f));
    % 
    % plot(time, this_peak_log_P_af_timecourse,'-k')
    % 
    % ylim([-10 20])
    % this_ylim=ylim;
    % hold on
    % plot([time(handles.ii_laser_start) time(handles.ii_laser_start)],[this_ylim(1) this_ylim(2)],'-k','LineWidth',2)
    % plot([time(handles.ii_laser_end) time(handles.ii_laser_end)],[this_ylim(1) this_ylim(2)],'-k','LineWidth',2)
    % 
    % title([handles.electrode_label{handles.peakLFPNo} ' Wavelet delta power timecourse for bandwidth from ' num2str(show_f_bandwidth(1)) ' to ' num2str(show_f_bandwidth(2)) ' Hz'])
    % xlabel('Time (min)')
    % ylabel('dB')
    % 

end

pffft=1
function handles_out=drgOBtoHIP_tPRP(handles_choices)
%drgOBChr2LFPpower 
%provides power analysis for Joe's ChR2 runs

if exist('handles_choices')==0
    clear all
    close all

    %     handles.peakLFPNo=9;%This is the LFP number that will be processed
    %     electrode_label='Right OB';

    %     handles.peakLFPNo=1;%This is the LFP number that will be processed
    %     electrode_label='Right CA1';

    handles.peakLFPNo=16;%This is the LFP number that will be processed
    handles.electrode_label{1}='Right HIP';
    handles.electrode_label{8}='Left HIP';
    handles.electrode_label{9}='Right OB';
    handles.electrode_label{16}='Left OB';
    electrode_label=handles.electrode_label{handles.peakLFPNo};

    %Choices for PAC frequency bandwidth
    handles.peakLowF=1;         %peak is the reference low frequency, usuallyy 6-14 Hz (theta) 
    handles.peakHighF=3;

    handles.burstLowF=30;       %burst is the power bandwidth (e.g. 65-95 Hz is high gamma)
    handles.burstHighF=50;

    % handles.peakLFPNo=1;%This is the LFP number that will be processed
    % electrode_label='Left CA1';

    prefix_save_file='prp_hgt'; %Always use prp_

    % handles.peakLFPNo=[1 8 9 16];
    show_f_bandwidth=[35 55];

  

    handles.n_phase_bins=50;

    handles.which_method=1;
    handles.save_drgb=0;
    handles.use_peakAngle=0;

    handles.window=1; %This is the FFT window in sec
    handles.noverlap=handles.window*0.9;

    % %Load file
    % handles.jtPathNames{1}='/Users/restrepd/Documents/Projects/Joe_OB_to_hippo/CNO/20231213_Siegfried_C5957_DREADD_n2_NewCNO_Pre_During_2mg_kg_2-Undecanone/20231213_Siegfried_C5957_DREADD_n2_NewCNO_Pre_231213_105533/';
    % handles.jtFileNames{1}='jt_times_20231213_Siegfried_C5957_DREADD_n2_NewCNO_Pre_231213_105533.mat';
    % handles.jtPathNames{2}='/Users/restrepd/Documents/Projects/Joe_OB_to_hippo/CNO/20231213_Siegfried_C5957_DREADD_n2_NewCNO_Pre_During_2mg_kg_2-Undecanone/20231213_Siegfried_C5957_DREADD_n2_NewCNO_2mg_kg_2-Undecanone_231213_111424/';
    % handles.jtFileNames{2}='jt_times_20231213_Siegfried_C5957_DREADD_n2_NewCNO_2mg_kg_2-Undecanone_231213_111424.mat';
    % 

    % handles.jtPathNames{1}='/Volumes/Diego HD/Joe/Optogenetics/5xFADvsWT_1_hour_treatmentVsnone/Curley_WT_Treated/20240206_WT_Curley_Tx_2/20240206_WT_Curley_Tx_2_240206_120854/';
    % handles.jtFileNames{1}='jt_times_20240206_WT_Curley_Tx_2_240206_120854.mat';

    % handles.jtPathNames{1}='/Volumes/Diego HD/Joe/Optogenetics/5xFADvsWT_1_hour_treatmentVsnone/Curley_WT_Treated/20240205_WT_Curley_Tx_1/20240205_WT_Curley_Tx_1_240205_101405/';
    % handles.jtFileNames{1}='jt_times_20240205_WT_Curley_Tx_1_240205_101405.mat'

% 
%     [jtFileName,jtPathName] = uigetfile('jt_times*.mat','Select jt_times file to open');
%     handles.jtfullName=[jtPathName,jtFileName];
    % handles.jtFileNames=jtFileNames;
    % handles.jtPathNames=jtPathNames;

    handles.showData=1;
    

    % ii_laser_start=60*5;
    % ii_laser_end=60*65;
    % 
    % handles.jtPathNames{1}='/Users/restrepd/Documents/Projects/Closed loop/Joe_ClosedLoop/20240306_Donatello_Closed_Loop_Test_withBandPass_4_240306_112305/';
    % handles.jtFileNames{1}='jt_times_20240306_Donatello_Closed_Loop_Test_withBandPass_4_240306_112305.mat';

     handles.jtPathNames{1}='/Users/restrepd/Documents/Projects/Joe_OB_to_hippo/5xFADvsWT_1_hour_treatmentVsnone/N_1/Curley_WT_Treated/20240205_WT_Curley_Tx_1/20240205_WT_Curley_Tx_1_240205_101405/';
    handles.jtFileNames{1}='jt_times_20240205_WT_Curley_Tx_1_240205_101405.mat';

    handles.ii_laser_start=34*9;
    handles.ii_laser_end=59*9;

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
    % handles.jtFileNames=jtFileNames;
    % handles.jtPathNames=jtPathNames;
  
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
 




handles.burstLFPNo=handles.peakLFPNo;

%Note: for contnuous
% t_start is 3.2 sec after the start of the trial and
% trials are spaced by 9 sec
%These times are set so that the entire trial is processed
%These settings read 9 sec
handles.time_start=-3.2-0.2;
handles.time_end=(9-3.2)+0.2;

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
handles_out.mod_indx=[];
handles_out.all_phase_histo=[];
handles_out.peakAngle=[];
handles_out.troughAngle=[];
handles_out.peakAngle_conv=[];
handles_out.troughAngle_conv=[];
handles_out.peakPowerPAC=[];
handles_out.troughPowerPAC=[];
handles_out.all_theta_wave=[];
time=[];
tNum_all=0;
ii_start_per_file=[];
ii_end_per_file=[];
for fileNo=1:length(handles.jtFileNames)

    jtPathName=handles.jtPathNames{fileNo};
    jtFileName=handles.jtFileNames{fileNo};
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



    [t_apt,freq,all_Power, all_Power_ref, all_Power_timecourse, this_trialNo, perCorr_pertr, which_event]=drgGetLFPwavePowerForThisEvTypeNo(handles);

    handles.drgb.PACwave.trialNos_PRP=this_trialNo;
    dt=handles.dt_tPRP;
    t_pac=[handles.time_start+handles.time_pad:dt:handles.time_end-handles.time_pad];

    %Get PAC
    % start_toc=toc;
    handles=drgThetaAmpPhaseTrialRange(handles);


    %Calculate a circ_mean peak angle
    no_conv_points=21;
    rad_meanPeakAngleConv=zeros(1,length(handles.drgb.PAC.peakAngle));

    rad_meanPeakAngle=pi*(handles.drgb.PAC.peakAngle-180)/180;
    %First half of the conv window
    rad_meanPeakAngleConv(1:((no_conv_points-1)/2))=circ_mean(rad_meanPeakAngle(1:((no_conv_points-1)/2))');

    %Last half of the conv window
    rad_meanPeakAngleConv(end-((no_conv_points-1)/2)+1:end)=circ_mean(rad_meanPeakAngle(end-((no_conv_points-1)/2)+1:end)');

    for ii=((no_conv_points-1)/2)+1:length(rad_meanPeakAngle)-((no_conv_points-1)/2)
        rad_meanPeakAngleConv(ii)=circ_mean(rad_meanPeakAngle(1,ii-((no_conv_points-1)/2):ii+((no_conv_points-1)/2))');
    end

    handles.drgb.PAC.meanPeakAngle_conv=(180/pi)*rad_meanPeakAngleConv+180;
 
    %Calculate a circ_mean trough angle
    rad_meanTroughAngleConv=zeros(1,length(handles.drgb.PAC.troughAngle));

    rad_meanTroughAngle=pi*(handles.drgb.PAC.troughAngle-180)/180;
    %First half of the conv window
    rad_meanTroughAngleConv(1:((no_conv_points-1)/2))=circ_mean(rad_meanTroughAngle(1:((no_conv_points-1)/2))');

    %Last half of the conv window
    rad_meanTroughAngleConv(end-((no_conv_points-1)/2)+1:end)=circ_mean(rad_meanTroughAngle(end-((no_conv_points-1)/2)+1:end)');

    for ii=((no_conv_points-1)/2)+1:length(rad_meanTroughAngle)-((no_conv_points-1)/2)
        rad_meanTroughAngleConv(ii)=circ_mean(rad_meanTroughAngle(1,ii-((no_conv_points-1)/2):ii+((no_conv_points-1)/2))');
    end

    handles.drgb.PAC.meanTroughAngle_conv=(180/pi)*rad_meanTroughAngleConv+180;

    handles_out.mod_indx=[handles_out.mod_indx handles.drgb.PAC.mod_indx];
    handles_out.all_phase_histo=[handles_out.all_phase_histo; handles.drgb.PAC.all_phase_histo];
    handles_out.out_times_env=handles.drgb.PAC.out_times_env;
    handles_out.phase=handles.drgb.PAC.phase;
    handles_out.peakAngle=[handles_out.peakAngle handles.drgb.PAC.peakAngle];
    handles_out.troughAngle=[handles_out.troughAngle handles.drgb.PAC.troughAngle];
    handles_out.peakAngle_conv=[handles_out.peakAngle_conv handles.drgb.PAC.meanPeakAngle_conv];
    handles_out.troughAngle_conv=[handles_out.troughAngle_conv handles.drgb.PAC.meanTroughAngle_conv];
    handles_out.peakPowerPAC=[handles_out.peakPowerPAC handles.drgb.PAC.meanPeakPower];
    handles_out.troughPowerPAC=[handles_out.troughPowerPAC handles.drgb.PAC.meanTroughPower];
    handles_out.all_theta_wave=[handles_out.all_theta_wave; handles.drgb.PAC.all_theta_wave];

    if fileNo==1
        ii_start_per_file(fileNo)=1;
    else
        ii_start_per_file(fileNo)=ii_end_per_file(fileNo-1)+1;
    end
    ii_end_per_file(fileNo)=length(handles_out.peakPowerPAC);

    % fprintf(1, 'drgThetaAmpPhaseTrialRange = %d\n',toc-start_toc)

    % start_toc=toc;
    handles.drgb.PACwave.trialNos_PAC=handles.drgb.PAC.this_trialNo;
    %If there was a problem with signal saturation for this LFP the number of
    %trials is zero, and we should skip further analysis

    if handles.drgb.PAC.no_trials>0
        %Please note that more trials are excluded from the PAC analysis than from
        %the lick analysis

        %Find the wavelet power at the peak and trough
        %     dt=handles.window-handles.noverlap;

        handles.drgb.PACwave.t_pac=t_pac;
        peakPower=zeros(handles.drgb.PAC.no_trials,length(t_pac));
        allPower=zeros(handles.drgb.PAC.no_trials,length(t_pac));
        %Find the value of gamma power at each point
        out_times_env=handles.drgb.PAC.out_times_env;
        out_times_env=out_times_env+handles.time_start+handles.time_pad;
        % meanPeakAngle=(handles.drgb.PAC.peakAngleForPower*(pi/180))-pi;
        % meanTroughAngle=(handles.drgb.PAC.troughAngleForPower*(pi/180))-pi;


        %Find peak wavelet power
        max_t_points=ceil(20*(t_apt(end)-t_apt(1)));
        % ii_start_per_file(fileNo)=tNum_all+1;

        %Find trough power
        troughPower=zeros(handles.drgb.PAC.no_trials,length(t_pac));

        for trNum=1:handles.drgb.PAC.no_trials

            tNum_all=tNum_all+1;


            if tNum_all==1
                time(1)=0;
            else
                time(tNum_all)=time(tNum_all-1)+handles.delta_trial;
            end

            if (trNum==1)&(fileNo==2)
                time(tNum_all)=time(tNum_all)+handles.delta_trial*4*60/9;
            end

            this_peakPower=zeros(1,max_t_points);
            this_peakPower_spectrum=zeros(max_t_points,length(freq));
            this_peakPower_times=zeros(1,max_t_points);

            ii=1;
            jj=0;
            at_end=0;
            this_angleTetaLFP=handles.drgb.PAC.PACtimecourse(trNum).decanglethetaLFP;
            %         this_LFPenv=handles.drgb.PAC.PACtimecourse(trNum).decLFPgenv;

            this_meanPeakAngle_conv=(handles.drgb.PAC.meanPeakAngle_conv(trNum)-180)*(pi/180);

            while at_end==0
                ii_next=find(this_angleTetaLFP(ii:end)>=this_meanPeakAngle_conv,1,'first');
                if (~isempty(ii_next))&(ii+ii_next-1<=length(t_apt))&(ii+ii_next-1<=length(out_times_env))
                    jj=jj+1;
                    this_peakPower(jj)=mean(10*log10(all_Power_timecourse(trNum,:,ii+ii_next-1)),2);
                    % this_peakPower_spectrum(jj,:)=10*log10(all_Power_timecourse(trNum,:,ii+ii_next-1));
                    this_peakPower_times(jj)=out_times_env(ii+ii_next-1);
                    ii=ii+ii_next;
                    ii_next=find(this_angleTetaLFP(ii:end)<this_meanPeakAngle_conv,1,'first');
                    if ~isempty(ii_next)
                        ii=ii+ii_next;
                    else
                        at_end=1;
                    end

                else
                    at_end=1;
                end
            end


            this_peakPower=this_peakPower(1,1:jj);
            % this_peakPower_spectrum=this_peakPower_spectrum(1:jj,:);
            this_peakPower_times=this_peakPower_times(1,1:jj);

            handles.drgb.PACwave.PACtimecourse(tNum_all).Power_per_peak=this_peakPower;
            handles.drgb.PACwave.PACtimecourse(tNum_all).peakPower_times=this_peakPower_times;
            handles_out.mean_peakPower(tNum_all)=mean(this_peakPower);


%             %         d_t_pac=t_pac(2)-t_pac(1);
% 
%             for ii_t=1:length(t_pac)
%                 if t_pac(ii_t)<=this_peakPower_times(1)
%                     peakPower(trNum,ii_t)=this_peakPower(1);
%                     %this_peakPower_t_pac(ii_t)=this_peakPower(1);
%                 else
%                     if t_pac(ii_t)>=this_peakPower_times(end)
%                         peakPower(trNum,ii_t)=this_peakPower(end);
%                         %this_peakPower_t_pac(ii_t)=this_peakPower(end);
%                     else
%                         ii_pt=find(this_peakPower_times>=t_pac(ii_t),1,'first');
%                         peakPower(trNum,ii_t)=this_peakPower(ii_pt-1)+(t_pac(ii_t)-this_peakPower_times(ii_pt-1))*((this_peakPower(ii_pt)-this_peakPower(ii_pt-1))/(this_peakPower_times(ii_pt)-this_peakPower_times(ii_pt-1)));
%                     end
%                 end
%             end
% 
% %             if handles.subtractRef==1
% %                 peakPower(trNum,:)=peakPower(trNum,:)-mean(peakPower(trNum,(t_pac>=handles.startRef+handles.time_pad)&(t_pac<=handles.endRef-handles.time_pad)));
% %             end
% 
%             handles.drgb.PACwave.peakPowerSpectrum(tNum_all,:)=mean(this_peakPower_spectrum,1);
%             handles.drgb.PACwave.PACtimecourse(tNum_all).peakPower=peakPower(trNum,:);
%             handles.drgb.PACwave.meanPeakPower(tNum_all)=mean(peakPower(trNum,:),2);

            %Calculate the power timecourse
            this_all_Power=zeros(1,length(t_pac));
            this_all_Power(1)=mean(mean(10*log10(all_Power_timecourse(trNum,:,t_apt<t_pac(1)+(dt/2))),2));
            this_all_Power(end)=mean(mean(10*log10(all_Power_timecourse(trNum,:,t_apt>t_pac(end)-(dt/2))),2));
            for ii=2:length(t_pac)-1
                this_all_Power(ii)=mean(mean(10*log10(all_Power_timecourse(trNum,:,(t_apt<t_pac(ii)+(dt/2))&(t_apt<t_pac(ii)+(dt/2)))),2));
            end

%             if handles.subtractRef==1
%                 this_all_Power=this_all_Power-mean(this_all_Power((t_pac>=handles.startRef+handles.time_pad)&(t_pac<=handles.endRef-handles.time_pad)));
%             end

            handles.drgb.PACwave.PACtimecourse(tNum_all).allPower=this_all_Power;
            allPower(trNum,:)=this_all_Power;

        % end
        % ii_end_per_file(fileNo)=tNum_all;
        % 
        % 
        % for trNum=1:handles.drgb.PAC.no_trials

            this_troughPower=zeros(1,max_t_points);
            % this_troughPower_spectrum=zeros(max_t_points,length(freq));
            this_troughPower_times=zeros(1,max_t_points);

            handles.drgb.PACwave.PACtimecourse(tNum_all).Power_per_trough=this_troughPower;
            handles.drgb.PACwave.PACtimecourse(tNum_all).troughPower_times=this_troughPower_times;

            ii=1;
            jj=0;
            at_end=0;
            this_angleTetaLFP=handles.drgb.PAC.PACtimecourse(trNum).decanglethetaLFP;
            %         this_LFPenv=handles.drgb.PAC.PACtimecourse(trNum).decLFPgenv;

            this_meanTroughAngle_conv=(handles.drgb.PAC.meanTroughAngle_conv(trNum)-180)*(pi/180);

            while at_end==0
                ii_next=find(this_angleTetaLFP(ii:end)>=this_meanTroughAngle_conv,1,'first');
                if (~isempty(ii_next))&(ii+ii_next-1<=length(t_apt))&(ii+ii_next-1<=length(out_times_env))
                    jj=jj+1;
                    this_troughPower(jj)=mean(10*log10(all_Power_timecourse(trNum,:,ii+ii_next-1)),2);
                    % this_troughPower_spectrum(jj,:)=10*log10(all_Power_timecourse(trNum,:,ii+ii_next-1));
                    this_troughPower_times(jj)=out_times_env(ii+ii_next-1);
                    ii=ii+ii_next;
                    ii_next=find(this_angleTetaLFP(ii:end)<this_meanTroughAngle_conv,1,'first');
                    if ~isempty(ii_next)
                        ii=ii+ii_next;
                    else
                        at_end=1;
                    end
                else
                    at_end=1;
                end
            end

            this_troughPower=this_troughPower(1,1:jj);
            % this_troughPower_spectrum=this_troughPower_spectrum(1:jj,:);
            this_troughPower_times=this_troughPower_times(1,1:jj);

            handles.drgb.PACwave.PACtimecourse(tNum_all).Power_per_trough=this_troughPower;
            handles.drgb.PACwave.PACtimecourse(tNum_all).troughPower_times=this_troughPower_times;
            handles_out.mean_troughPower(tNum_all)=mean(this_troughPower);

            % for ii_t=1:length(t_pac)
            %     if t_pac(ii_t)<=this_troughPower_times(1)
            %         troughPower(trNum,ii_t)=this_troughPower(1);
            %         %this_troughPower_t_pac(ii_t)=this_troughPower(1);
            %     else
            %         if t_pac(ii_t)>=this_troughPower_times(end)
            %             troughPower(trNum,ii_t)=this_troughPower(end);
            %             %this_troughPower_t_pac(ii_t)=this_troughPower(end);
            %         else
            %             ii_pt=find(this_troughPower_times>=t_pac(ii_t),1,'first');
            %             troughPower(trNum,ii_t)=this_troughPower(ii_pt-1)+(t_pac(ii_t)-this_troughPower_times(ii_pt-1))*((this_troughPower(ii_pt)-this_troughPower(ii_pt-1))/(this_troughPower_times(ii_pt)-this_troughPower_times(ii_pt-1)));
            %         end
            %     end
            % end

%             if handles.subtractRef==1
%                 troughPower(trNum,:)=troughPower(trNum,:)-mean(troughPower(trNum,(t_pac>=handles.startRef+handles.time_pad)&(t_pac<=handles.endRef-handles.time_pad)));
%             end

            % handles.drgb.PACwave.troughPowerSpectrum(tNum_all,:)=mean(this_troughPower_spectrum,1);
            % handles.drgb.PACwave.PACtimecourse(tNum_all).troughPower=troughPower(trNum,:);
            % handles.drgb.PACwave.meanTroughPower(tNum_all)=mean(troughPower(trNum,:),2);
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
end
%
% ii_laser_start=floor(find(dec_laser>0.5,1,'first')/(1000*mean_window));
% ii_laser_end=ceil(find(dec_laser>0.5,1,'last')/(1000*mean_window));

handles_out.time=time;
handles_out.ii_start_per_file=ii_start_per_file;
handles_out.ii_end_per_file=ii_end_per_file;

save_base=handles.jtFileNames{1};
save([handles.jtPathNames{1} prefix_save_file '_ch' num2str(handles.peakLFPNo), '_', save_base(9:end)],'handles','handles_out')

% close all
figNo=8;

time=time/60; %Change time to minutes

%Plot Modulation index
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig2 = figure(figNo);
set(hFig2, 'units','normalized','position',[.25 .1 .4 .25])
hold on

mod_indx=handles_out.mod_indx;

plot(time,mod_indx,'ob-')


if length(mod_indx)>20
    no_conv_points=6;
    conv_win=ones(1,no_conv_points);
    mod_indx_extend=[mean(mod_indx(1:no_conv_points/2))*ones(1,no_conv_points) mod_indx mean(mod_indx(end-(no_conv_points/2)+1:end))*ones(1,no_conv_points)];
    conv_mod_indx_extend=conv(mod_indx_extend,conv_win)/no_conv_points;
    conv_mod_indx=conv_mod_indx_extend(no_conv_points+1:no_conv_points+length(mod_indx));
    plot(time,conv_mod_indx,'-b','LineWidth',3)
end
theseyl=ylim;

xlabel('Time (min)')
ylabel('Modulation index')
title('Modulation index vs time')

%Peak and trough processed above
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig2 = figure(figNo);
set(hFig2, 'units','normalized','position',[.25 .1 .4 .25])
hold on

reference_power=(   mean(handles_out.mean_troughPower(time<handles.ii_laser_start/60))  +...
    mean(handles_out.mean_peakPower(time<handles.ii_laser_start/60)) )/2;
meanTroughPower=handles_out.mean_troughPower-reference_power;
meanPeakPower=handles_out.mean_peakPower-reference_power;

plot(time(ii_start_per_file(1):ii_end_per_file(1)),meanTroughPower(ii_start_per_file(1):ii_end_per_file(1)),'b-')
% plot(time(ii_start_per_file(2):ii_end_per_file(2)),meanTroughPower(ii_start_per_file(2):ii_end_per_file(2)),'b-')

plot(time(ii_start_per_file(1):ii_end_per_file(1)),meanPeakPower(ii_start_per_file(1):ii_end_per_file(1)),'r-')
% plot(time(ii_start_per_file(2):ii_end_per_file(2)),meanPeakPower(ii_start_per_file(2):ii_end_per_file(2)),'r-')

% plot([between_ii between_ii],theseyl,'-k','LineWidth',3)

this_xl=xlim;
plot([this_xl],[0 0],'-k')

xlabel('Time (min)')
ylabel('Peak power (dB)')


legend('Through','Peak')
title(['Power (dB)  '])

%Peak and trough processed through PAC
%Comenting this because I do not see differences with the graph above
% figNo=figNo+1;
% try
%     close(figNo)
% catch
% end
% 
% hFig2 = figure(figNo);
% set(hFig2, 'units','normalized','position',[.25 .1 .4 .25])
% hold on
% 
% reference_power=(mean(handles_out.troughPowerPAC(ii_start_per_file(1):ii_end_per_file(1)))+...
%     mean(handles_out.peakPowerPAC(ii_start_per_file(1):ii_end_per_file(1))))/2;
% meanTroughPower=handles_out.troughPowerPAC-reference_power;
% meanPeakPower=handles_out.peakPowerPAC-reference_power;
% 
% plot(time(ii_start_per_file(1):ii_end_per_file(1)),meanTroughPower(ii_start_per_file(1):ii_end_per_file(1)),'b-')
% % plot(time(ii_start_per_file(2):ii_end_per_file(2)),meanTroughPower(ii_start_per_file(2):ii_end_per_file(2)),'b-')
% 
% plot(time(ii_start_per_file(1):ii_end_per_file(1)),meanPeakPower(ii_start_per_file(1):ii_end_per_file(1)),'r-')
% % plot(time(ii_start_per_file(2):ii_end_per_file(2)),meanPeakPower(ii_start_per_file(2):ii_end_per_file(2)),'r-')
% 
% % plot([between_ii between_ii],theseyl,'-k','LineWidth',3)
% 
% this_xl=xlim;
% plot([this_xl],[0 0],'-k')
% 
% xlabel('Time(min)')
% ylabel('Peak power (dB)')
% 
% 
% legend('Through','Peak')
% title(['Power (dB, PAC calculated)  '])


%Peak angle
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig2 = figure(figNo);
set(hFig2, 'units','normalized','position',[.25 .2 .7 .23])
hold on


peakAngle=handles_out.peakAngle;
peakAngle_conv=handles_out.peakAngle_conv;

troughAngle=handles_out.troughAngle;
troughAngle_conv=handles_out.troughAngle_conv;

plot(time(ii_start_per_file(1):ii_end_per_file(1)),troughAngle(ii_start_per_file(1):ii_end_per_file(1)),'bo-')
% plot(time(ii_start_per_file(2):ii_end_per_file(2)),troughAngle(ii_start_per_file(2):ii_end_per_file(2)),'bo-')

plot(time(ii_start_per_file(1):ii_end_per_file(1)),troughAngle_conv(ii_start_per_file(1):ii_end_per_file(1)),'b-','LineWidth',3)
% plot(time(ii_start_per_file(2):ii_end_per_file(2)),troughAngle_conv(ii_start_per_file(2):ii_end_per_file(2)),'b-','LineWidth',3)

plot(time(ii_start_per_file(1):ii_end_per_file(1)),peakAngle(ii_start_per_file(1):ii_end_per_file(1)),'ro-')
% plot(time(ii_start_per_file(2):ii_end_per_file(2)),peakAngle(ii_start_per_file(2):ii_end_per_file(2)),'ro-')

plot(time(ii_start_per_file(1):ii_end_per_file(1)),peakAngle_conv(ii_start_per_file(1):ii_end_per_file(1)),'r-','LineWidth',3)
% plot(time(ii_start_per_file(2):ii_end_per_file(2)),peakAngle_conv(ii_start_per_file(2):ii_end_per_file(2)),'r-','LineWidth',3)

% plot([between_ii between_ii],theseyl,'-k','LineWidth',3)

xlabel('Time(min)')
ylabel('Phase')


legend('Through','Peak')
title(['Angle   '])


figNo=figNo+1;
try
    close(figNo)
catch
end

hFig2 = figure(figNo);
set(hFig2, 'units','normalized','position',[.25 .2 .7 .23])
hold on

all_phase_histo=handles_out.all_phase_histo;

min_prob=prctile(all_phase_histo(:),5);
max_prob=prctile(all_phase_histo(:),95);


%pcolor(repmat(phase,length(trials),1),repmat(trials',1,length(phase)),all_phase_histo)
drg_pcolor(repmat(time,length(handles_out.phase),1),repmat(handles_out.phase',1,length(time)),all_phase_histo')
%         colormap jet
colormap fire
shading flat
% min_prob=0.0113;
% max_prob=0.0314;
caxis([min_prob    max_prob])
xlabel('Time (min)')
ylabel('Phase for low freq oscillation (deg)');
title(['Phase-amplitude coupling timecourse '])

figNo=figNo+1;
try
    close(figNo)
catch
end

hFig2 = figure(figNo);
set(hFig2, 'units','normalized','position',[.25 .25 .25 .25])
hold on

shadedErrorBar(handles.drgb.PAC.phase,mean(handles_out.all_theta_wave,1),std(handles_out.all_theta_wave,0,1),'-b')
xlim([0 360])
title('Mean low frequency waveform')
xlabel('Degrees')

pffft=1
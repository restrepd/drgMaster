function handles=drgLFPwaveSpectrogram(handles)
%Generates a timecourse of the LFP power in decibels 10*log10(Power)

% tic

% first_toc=toc;

handles.drgb.PACwave=[];
handles.drgb.lickStuff=[];

if isfield(handles,'calculate_lick')
    calculate_lick=handles.calculate_lick;
else
    %The default is to calculate the licks
    calculate_lick=1;
end
 
% start_toc=toc;
[t_apt,freq,all_Power,all_Power_ref, all_Power_timecourse, this_trialNo, perCorr_pertr, which_event]=drgGetLFPwavePowerForThisEvTypeNo(handles);
% fprintf(1, 'dt for drgGetLFPwavePowerForThisEvTypeNo = %d\n',toc-start_toc)

if handles.displayData==1
    if ~isempty(this_trialNo)

        %Timecourse doing average after log
        %Get max and min
        if handles.subtractRef==0
            log_P_timecourse=zeros(length(freq),length(t_apt));
            log_P_timecourse(:,:)=mean(10*log10(all_Power_timecourse),1);

            %Per trial power plot
            log_P_per_trial_timecourse=zeros(length(freq)*length(this_trialNo),length(t_apt));
            y_shift=0;
            for trialNo=1:length(this_trialNo)
                this_log_P_timecourse=zeros(length(freq),length(t_apt));
                this_log_P_timecourse(:,:)=10*log10(all_Power_timecourse(trialNo,:,:));
                log_P_per_trial_timecourse(y_shift+1:y_shift+length(freq),:)=this_log_P_timecourse;
                shifted_freq(1,y_shift+1:y_shift+length(freq))=freq+(trialNo-1)*freq(end);
                y_shift=y_shift+length(freq);
            end

            if handles.autoscale==1
                maxLogPper=prctile(log_P_timecourse(:),99);
                minLogPper=prctile(log_P_timecourse(:),1);
                %Note: Diego added this on purpose to limit the range to 10 dB
                %This results in emphasizing changes in the top 10 dB
                if maxLogPper-minLogPper>12
                    minLogPper=maxLogPper-12;
                end
            else
                maxLogPper=handles.maxLogP;
                minLogPper=handles.minLogP;
            end
        else
            log_P_timecourse=zeros(length(freq),length(t_apt));
            log_P_timecourse(:,:)=mean(10*log10(all_Power_timecourse),1);
            log_P_timecourse_ref=zeros(length(freq),length(t_apt));
            log_P_timecourse_ref(:,:)=repmat(mean(10*log10(all_Power_ref),1)',1,length(t_apt));

            %Per trial power plot
            log_P_per_trial_timecourse_sub=zeros(length(freq)*length(this_trialNo),length(t_apt));
            y_shift=0;
            sy_shift=0;
            shifted_freq=[];
            for trialNo=1:length(this_trialNo)
                this_log_P_timecourse=zeros(length(freq),length(t_apt));
                this_log_P_timecourse(:,:)=10*log10(all_Power_timecourse(trialNo,:,:));
                this_log_P_timecourse_ref=zeros(length(freq),length(t_apt));
                this_log_P_timecourse_ref(:,:)=repmat(mean(10*log10(all_Power_ref(trialNo,:)),1)',1,length(t_apt));
                log_P_per_trial_timecourse_sub(y_shift+1:y_shift+length(freq),:)=this_log_P_timecourse-this_log_P_timecourse_ref;
                shifted_freq(1,y_shift+1:y_shift+length(freq))=freq+(trialNo-1)*freq(end);
                y_shift=y_shift+length(freq);
            end

            max_delta=16;
            if handles.autoscale==1

                deltaLogP=log_P_timecourse'-log_P_timecourse_ref';
                maxLogPper=prctile(deltaLogP(:),99);
                minLogPper=prctile(deltaLogP(:),1);
                %Note: Diego added this on purpose to limit the range to 10 dB
                %This results in emphasizing changes in the top 10 dB
                if maxLogPper-minLogPper>max_delta
                    minLogPper=maxLogPper-max_delta;
                end

            else
                maxLogPper=handles.maxLogP;
                minLogPper=handles.minLogP;
            end
        end

        figNo=0;

        figNo=figNo+1;
        try
            close(figNo)
        catch
        end

        %Plot the timecourse
        hFig = figure(figNo);
        set(hFig, 'units','normalized','position',[.07 .05 .75 .3])
        if handles.subtractRef==0
            drg_pcolor(repmat(t_apt,length(freq),1)',repmat(freq,length(t_apt),1),log_P_timecourse')
        else
            %pcolor(repmat(t,length(f),1)',repmat(f,length(t),1),10*log10(P_timecourse')-10*log10(P_timecourse_ref'))
            drg_pcolor(repmat(t_apt,length(freq),1)',repmat(freq,length(t_apt),1),log_P_timecourse'-log_P_timecourse_ref')
            %imagesc(t,f,10*log10(P_timecourse')-10*log10(P_timecourse_ref'))
        end

        colormap fire
        shading interp
        caxis([minLogPper maxLogPper]);
        xlabel('Time (sec)')
        ylabel('Frequency (Hz)');
        title(['Power (dB, wavelet) timecourse ' handles.drg.session(1).draq_d.eventlabels{handles.evTypeNo}])


        figNo=figNo+1;
        try
            close(figNo)
        catch
        end

        hFig = figure(figNo);
        set(hFig, 'units','normalized','position',[.83 .1 .05 .3])

        prain=[minLogPper:(maxLogPper-minLogPper)/99:maxLogPper];
        drg_pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
        colormap fire
        shading interp
        ax=gca;
        set(ax,'XTickLabel','')

        %Plot the dependence on frequency
        figNo=figNo+1;
        hFig = figure(figNo);
        set(hFig, 'units','normalized','position',[.1 .1 .4 .4])
        hold on


        this_mean_dbWB=zeros(1,length(freq));
        this_mean_dbWB(1,:)=mean(log_P_timecourse',1);

        try
            CI=[];
            CI = bootci(1000, {@mean, log_P_timecourse'})';
            CI(:,1)= this_mean_dbWB'-CI(:,1);
            CI(:,2)=CI(:,2)- this_mean_dbWB';


            [hlCR, hpCR] = boundedline(freq',this_mean_dbWB', CI, 'b');
        catch
            plot(freq',this_mean_dbWB', 'b');
        end


        title('Wavelet power spectrum')
        xlabel('Frequency (Hz)')
        ylabel('dB')


    end
end

pfft=1;



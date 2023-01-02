function drgDisplayBatchMultiDayLFP

close all
clear all

[FileName,PathName] = uigetfile({'*.mat'},'Select the .mat file with power analysis');
fprintf(1, ['\ndrgDisplayBatchMultiDayLFP run for ' FileName '\n\n']);

%Load the data
load([PathName FileName])

f_label{1}='Delta';
f_label{2}='Theta';
f_label{3}='Beta';
f_label{4}='Gamma_low';
f_label{5}='Gamma_high';

use_einf=0;

einfSignalLabels=[];
if use_einf==1
    einf=edfinfo([PathName handlespf.drg.drta_p.FileName(1:end-3) 'edf']);
    for chNo=1:4
        einfSignalLabels{chNo}=einf.SignalLabels{chNo};
    end
else
    einfSignalLabels{1}='Hippocampus left';
    einfSignalLabels{2}='Olfactory bulb';
    einfSignalLabels{3}='Cortex';
    einfSignalLabels{1}='Hippocampus right';
end

%Plot the data
figNo=0;



t_hours=handlespf.drgbchoices.start_hour+[1:size(mean_dB_power,2)]*(9/(60*60));

%Plot dB per bandwidth
for chNo=1:4
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);
    set(hFig, 'units','normalized','position',[.07 .05 .45 .75])
    sgtitle(['Power for ' einfSignalLabels{chNo}])
    for ii_plot=1:5
        subplot(5,1,ii_plot)
        this_mean_dB_power=zeros(1,size(mean_dB_power,2));
        this_mean_dB_power(1,:)=mean_dB_power(ii_plot,:,chNo);
        plot(t_hours,this_mean_dB_power)

        min_dB=prctile(this_mean_dB_power,1)-0.1*(prctile(this_mean_dB_power,99)-prctile(this_mean_dB_power,1));
        max_dB=prctile(this_mean_dB_power,99)+0.1*(prctile(this_mean_dB_power,99)-prctile(this_mean_dB_power,1));
        ylim([min_dB max_dB])
        title(f_label{ii_plot})
        ylabel('dB')
        xlabel('time (hrs)')


        t=handlespf.drgbchoices.light_on;

        %This only works if the start is during the day
        hold on
        while t< t_hours(end)

            if (rem(t,24)<handlespf.drgbchoices.light_on)||(rem(t,24)>=handlespf.drgbchoices.light_off)
                %This is the night

                delta_t_night_end=24-handlespf.drgbchoices.light_off+handlespf.drgbchoices.light_on;
                if t+delta_t_night_end<t_hours(end)
                    plot([t t+delta_t_night_end], [min_dB+0.1*(prctile(this_mean_dB_power,99)-prctile(this_mean_dB_power,1))...
                        min_dB+0.1*(prctile(this_mean_dB_power,99)-prctile(this_mean_dB_power,1)) ], '-k', 'LineWidth', 3)
                else
                    plot([t t_hours(end)], [min_dB+0.1*(prctile(this_mean_dB_power,99)-prctile(this_mean_dB_power,1))...
                        min_dB+0.1*(prctile(this_mean_dB_power,99)-prctile(this_mean_dB_power,1)) ], '-k', 'LineWidth', 3)
                end
                t=t+delta_t_night_end;

            else
                %This is the day
                t_last=t;
                t=t+(handlespf.drgbchoices.light_off-rem(t,24));


            end

        end

    end
end

%Subtract delta
%Plot dB per bandwidth
for chNo=1:4
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);
    set(hFig, 'units','normalized','position',[.07 .05 .45 .75])
    sgtitle(['Delta power subtracted for ' einfSignalLabels{chNo}])
    for ii_plot=2:5
        subplot(4,1,ii_plot-1)
        this_mean_dB_power=zeros(1,size(mean_dB_power,2));
        this_mean_dB_power(1,:)=mean_dB_power(ii_plot,:,chNo)-mean_dB_power(1,:,chNo);
        plot(t_hours,this_mean_dB_power)

        min_dB=prctile(this_mean_dB_power,1)-0.1*(prctile(this_mean_dB_power,99)-prctile(this_mean_dB_power,1));
        max_dB=prctile(this_mean_dB_power,99)+0.1*(prctile(this_mean_dB_power,99)-prctile(this_mean_dB_power,1));
        ylim([min_dB max_dB])
        title(f_label{ii_plot})
        ylabel('delta dB')
        xlabel('time (hrs)')

        t=handlespf.drgbchoices.light_on;

        %This only works if the start is during the day
        hold on
        while t< t_hours(end)

            if (rem(t,24)<handlespf.drgbchoices.light_on)||(rem(t,24)>=handlespf.drgbchoices.light_off)
                %This is the night

                delta_t_night_end=24-handlespf.drgbchoices.light_off+handlespf.drgbchoices.light_on;
                if t+delta_t_night_end<t_hours(end)
                    plot([t t+delta_t_night_end], [min_dB+0.1*(prctile(this_mean_dB_power,99)-prctile(this_mean_dB_power,1))...
                        min_dB+0.1*(prctile(this_mean_dB_power,99)-prctile(this_mean_dB_power,1)) ], '-k', 'LineWidth', 3)
                else
                    plot([t t_hours(end)], [min_dB+0.1*(prctile(this_mean_dB_power,99)-prctile(this_mean_dB_power,1))...
                        min_dB+0.1*(prctile(this_mean_dB_power,99)-prctile(this_mean_dB_power,1)) ], '-k', 'LineWidth', 3)
                end
                t=t+delta_t_night_end;

            else
                %This is the day
                t_last=t;
                t=t+(handlespf.drgbchoices.light_off-rem(t,24));


            end

        end

    end
end

pfft=1;
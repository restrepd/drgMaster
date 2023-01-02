%drgSummaryMultiDayLFPv2
close all
clear all

%Variables use to process the data
this_per=1; %Perecentile for maximum and minimum ylim
skip_ii=250; %Skip the first and last 10 sec because of large LFP changes 
figNo=0; %Initialize index for figure numbers

conv_window=500; %convolution window to smooth REM/NREM changes during sleep, I use 500
v=ones(1,conv_window)/conv_window;

%Below enter information for each file
PathName='/Users/restrepd/Documents/Projects/Multi day LFP/';
 
group_name{1}='Vehicle';
group_name{2}='VZV treated';
group_name{3}='Mock';

einfSignalLabels{1}='Hippocampus left';
einfSignalLabels{2}='Olfactory bulb';
einfSignalLabels{3}='Cortex';
einfSignalLabels{4}='Hippocampus right';

f_label{1}='Delta';
f_label{2}='Theta';
f_label{3}='Beta';
f_label{4}='Gamma low';
f_label{5}='Gamma high';

FileNames.mouse(1).FileName{1}='Restrepo C031_2022-10-26_10_00_14_export_ov0p5_fft.mat';
FileNames.mouse(1).date{1}='10/26/2022';
FileNames.mouse(1).FileName{2}='Restrepo C031_2_2022-11-03_11_22_23_export_ov0p5_fft.mat';
FileNames.mouse(1).date{2}='11/03/2022';
group(1)=1;
treatment_start(1)=152;

FileNames.mouse(2).FileName{1}='Restrepo T034_2022-10-26_10_00_29_export_ov0p5_fft.mat';
FileNames.mouse(2).date{1}='10/26/2022';
FileNames.mouse(2).FileName{2}='Restrepo T034_2_2022-11-03_11_35_29_export_ov0p5_fft.mat';
FileNames.mouse(2).date{2}='11/03/2022';
group(2)=2;
treatment_start(2)=152;

FileNames.mouse(3).FileName{1}='Restrepo T035_2022-10-26_10_00_19_export_ov0p5_fft.mat';
FileNames.mouse(3).date{1}='10/26/2022';
FileNames.mouse(3).FileName{2}='Restrepo T035_2_2022-11-03_11_27_21_export_ov0p5_fft.mat';
FileNames.mouse(3).date{2}='11/03/2022';
group(3)=3;
treatment_start(3)=152;

%Now process the data
LFPs=[];

%Plot dB per bandwidth
%The dB timecourse will be convolved with a flat window
for mouseNo=1:length(FileNames.mouse)
 
    t_hours=[];
    all_mean_dB_power=[];
    for fileNo=1:length(FileNames.mouse(mouseNo).FileName)

        load([PathName FileNames.mouse(mouseNo).FileName{fileNo}])
  
        numdays = datenum(FileNames.mouse(mouseNo).date{fileNo}) - datenum(FileNames.mouse(mouseNo).date{1});
        t_hours=[t_hours handlespf.drgbchoices.start_hour+24*numdays+[1:size(mean_dB_power,2)]*(9/(60*60))];
        all_mean_dB_power=[all_mean_dB_power mean_dB_power];
        LFP.mice(mouseNo).t_hours=t_hours;
        if fileNo==1
            t=handlespf.drgbchoices.start_hour;
        end
    end

    %This only works if the start is during the day
    LFP.mice(mouseNo).t_night_on=[];
    LFP.mice(mouseNo).t_night_off=[];
 
    ii_night=0;
    while t< t_hours(end)

        if ii_night==8
            pffft=1;
        end
        if (rem(t,24)<handlespf.drgbchoices.light_on)||(rem(t,24)>=handlespf.drgbchoices.light_off)

            %This is the night
            delta_t_night_end=24-handlespf.drgbchoices.light_off+handlespf.drgbchoices.light_on;
            ii_night=ii_night+1;
            if t+delta_t_night_end<t_hours(end)
                LFP.mice(mouseNo).t_night_on(ii_night)=t;
                LFP.mice(mouseNo).t_night_off(ii_night)=t+delta_t_night_end;
            else
                LFP.mice(mouseNo).t_night_on(ii_night)=t;
                LFP.mice(mouseNo).t_night_off(ii_night)=t_hours(end);
            end

            t=t+delta_t_night_end;

        else
            %This is the day
            t=t+(handlespf.drgbchoices.light_off-rem(t,24));
        end

    end

    for chNo=1:4

        for ii_plot=1:5
            this_mean_dB_power=zeros(1,size(all_mean_dB_power,2));
            this_mean_dB_power(1,:)=all_mean_dB_power(ii_plot,:,chNo);
            LFP.ch(chNo).bandwidth(ii_plot).mouse(mouseNo).mean_dB_power=this_mean_dB_power;
        end
    end


end






%Plot dB per bandwidth
%The dB timecourse will be convolved with a Gaussian
for chNo=1:4

    for ii_plot=1:5

        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
   
        hFig = figure(figNo);
        set(hFig, 'units','normalized','position',[.07 .05 .7 .7])
        hold on

        sgtitle(['Power for ' einfSignalLabels{chNo} ' ' f_label{ii_plot}])

        for grNo=1:max(group)

            subplot(max(group),1,grNo)
            hold on

            for mouseNo=1:length(LFP.ch(chNo).bandwidth(ii_plot).mouse)
                if group(mouseNo)==grNo
                    this_mean_dB_power=zeros(1,length(LFP.ch(chNo).bandwidth(ii_plot).mouse(mouseNo).mean_dB_power));
                    this_mean_dB_power(1,:)=LFP.ch(chNo).bandwidth(ii_plot).mouse(mouseNo).mean_dB_power;

                    conv_this_mean_dB_power=conv(this_mean_dB_power,v,'same');

                    %Zero dB by subtracting the average of the day before
                    t_zero_start=LFP.mice(mouseNo).t_night_off(floor(treatment_start(mouseNo)/24)-1);
                    t_zero_end=LFP.mice(mouseNo).t_night_on(floor(treatment_start(mouseNo)/24));
                    dB_zero=mean(conv_this_mean_dB_power((t_hours>=t_zero_start)&(t_hours<=t_zero_end)));
                    conv_this_mean_dB_power=conv_this_mean_dB_power-dB_zero;
                    skip_conv_this_mean_dB_power=conv_this_mean_dB_power(skip_ii+1:end-skip_ii);
                    min_dB=prctile(skip_conv_this_mean_dB_power,this_per)-0.1*(prctile(skip_conv_this_mean_dB_power,100-this_per)-prctile(skip_conv_this_mean_dB_power,this_per));
                    max_dB=prctile(skip_conv_this_mean_dB_power,100-this_per)+0.1*(prctile(skip_conv_this_mean_dB_power,100-this_per)-prctile(skip_conv_this_mean_dB_power,this_per));

                    %Plot dB power
                    t_hours=LFP.mice(mouseNo).t_hours;
                    plot(t_hours(skip_ii+1:end-skip_ii),conv_this_mean_dB_power(skip_ii+1:end-skip_ii),'-k')
                    
                    %Plot markers for night
                    for ii_nights=1:length(LFP.mice(mouseNo).t_night_on)
                        plot([LFP.mice(mouseNo).t_night_on(ii_nights) LFP.mice(mouseNo).t_night_off(ii_nights)],...
                            [min_dB+0.1*(prctile(this_mean_dB_power,99)-prctile(this_mean_dB_power,1))...
                            min_dB+0.1*(prctile(this_mean_dB_power,99)-prctile(this_mean_dB_power,1)) ], '-k', 'LineWidth', 3)

                        plot([LFP.mice(mouseNo).t_night_on(ii_nights) LFP.mice(mouseNo).t_night_on(ii_nights)],...
                            [min_dB max_dB], '-k')

                         plot([LFP.mice(mouseNo).t_night_off(ii_nights) LFP.mice(mouseNo).t_night_off(ii_nights)],...
                            [min_dB max_dB], '-k')

                    end
                    
                    %Place a vertical bar at start of treatment
                    plot([treatment_start(mouseNo) treatment_start(mouseNo)], [min_dB max_dB],'-r','LineWidth',3)


                end

            end


         
            title(group_name{grNo})

            ylim([min_dB max_dB])
            ylabel('dB')
            xlabel('time (hrs)')
        end


    end
end

fprintf(1, 'The LFP data were convolved within a window of %d hours\n',conv_window*(t_hours(2)-t_hours(1)));

pffft=1;
%drgSummaryMultiDayLFPv2
close all
clear all

%Variables use to process the data
this_per=1; %Perecentile for maximum and minimum ylim
skip_ii=250; %Skip the first and last 10 sec because of large LFP changes 
figNo=0; %Initialize index for figure numbers

conv_window=1; %convolution window to smooth REM/NREM changes during sleep, I use 500
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

sec_per_LFP_point=9;

fig_xlim=[80 240];
fig_ylim=[-3.5 5];
override_lims=0;

is_c57=1;

if is_c57==1

    %C57
    FileNames.mouse(1).FileName{1}='Restrepo C031_2022-10-26_10_00_14_export_ov0p5_fft.mat';
    FileNames.mouse(1).date{1}='10/26/2022';
    FileNames.mouse(1).FileName{2}='Restrepo C031_2_2022-11-03_11_22_23_export_ov0p5_fft.mat';
    FileNames.mouse(1).date{2}='11/03/2022';
    group(1)=1;
    FileNames.mouse(1).treatment_start=[152 152+24 152+48];

    FileNames.mouse(2).FileName{1}='Restrepo T034_2022-10-26_10_00_29_export_ov0p5_fft.mat';
    FileNames.mouse(2).date{1}='10/26/2022';
    FileNames.mouse(2).FileName{2}='Restrepo T034_2_2022-11-03_11_35_29_export_ov0p5_fft.mat';
    FileNames.mouse(2).date{2}='11/03/2022';
    group(2)=2;
    FileNames.mouse(2).treatment_start=[152 152+24 152+48];

    FileNames.mouse(3).FileName{1}='Restrepo T035_2022-10-26_10_00_19_export_ov0p5_fft.mat';
    FileNames.mouse(3).date{1}='10/26/2022';
    FileNames.mouse(3).FileName{2}='Restrepo T035_2_2022-11-03_11_27_21_export_ov0p5_fft.mat';
    FileNames.mouse(3).date{2}='11/03/2022';
    group(3)=3;
    FileNames.mouse(3).treatment_start=[152 152+24 152+48];

else

    %5xFAD
    FileNames.mouse(1).FileName{1}='Restrepo T022_2022-11-07_09_18_56_export_ov0p5_fft.mat';
    FileNames.mouse(1).date{1}='11/07/2022';
    FileNames.mouse(1).FileName{2}='Restrepo T022_2_2022-11-14_09_30_26_export_ov0p5_fft.mat';
    FileNames.mouse(1).date{2}='11/14/2022';
    FileNames.mouse(1).FileName{3}='Restrepo T022_3_2022-11-21_08_24_12_export_ov0p5_fft.mat';
    FileNames.mouse(1).date{3}='11/21/2022';
    FileNames.mouse(1).FileName{4}='Restrepo T022_4_2022-11-28_08_23_45_export_ov0p5_fft.mat';
    FileNames.mouse(1).date{4}='11/28/2022';
    group(1)=3;
    FileNames.mouse(1).treatment_start=[104.9 189.6 329.6];

    FileNames.mouse(2).FileName{1}='T032 22-11-07 export_ov0p5_fft.mat';
    FileNames.mouse(2).date{1}='11/07/2022';
    FileNames.mouse(2).FileName{2}='Restrepo T032_2_2022-11-14_09_30_19_export_ov0p5_fft.mat';
    FileNames.mouse(2).date{2}='11/14/2022';
    FileNames.mouse(2).FileName{3}='Restrepo T032_3_2022-11-21_08_24_04_export_ov0p5_fft.mat';
    FileNames.mouse(2).date{3}='11/21/2022';
    FileNames.mouse(2).FileName{4}='Restrepo T032_4_2022-11-28_08_23_38_export_ov0p5_fft.mat';
    FileNames.mouse(2).date{4}='11/28/2022';
    group(2)=2;
    FileNames.mouse(2).treatment_start=[104.9 189.6 329.6];

end

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
        start_hour=handlespf.drgbchoices.start_hour
        t_hours=[t_hours handlespf.drgbchoices.start_hour(1)+24*numdays+[1:size(mean_dB_power,2)]*sec_per_LFP_point/(60*60)];
        all_mean_dB_power=[all_mean_dB_power mean_dB_power];
        LFP.mice(mouseNo).t_hours=t_hours;
        if fileNo==1
            t=handlespf.drgbchoices.start_hour(1);
        end
    end

    %This only works if the start is during the day
    LFP.mice(mouseNo).t_night_on=[];
    LFP.mice(mouseNo).t_night_off=[];
 
    ii_night=0;
    ii_delta=0;
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

            %Do the delta night - day calculations
            if ii_night>=2
                ii_delta=ii_delta+1;
                for chNo=1:4
                    for ii_plot=1:5
                        this_mean_dB_power=zeros(1,size(all_mean_dB_power,2));
                        this_mean_dB_power(1,:)=all_mean_dB_power(ii_plot,:,chNo);
                        t_night_start=LFP.mice(mouseNo).t_night_on(ii_night);
                        t_night_end=LFP.mice(mouseNo).t_night_off(ii_night);
                        t_day_start=LFP.mice(mouseNo).t_night_off(ii_night-1);
                        t_day_end=LFP.mice(mouseNo).t_night_on(ii_night);
                        LFP.ch(chNo).bandwidth(ii_plot).mouse(mouseNo).deltadBtimepoint(ii_delta).delta_dB_power_night_day=mean(this_mean_dB_power((t_hours>=t_night_start)&(t_hours<=t_night_end)))-...
                            mean(this_mean_dB_power((t_hours>=t_day_start)&(t_hours<=t_day_end)));
                    end
                end
            end

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

            mouse_here=0;
            for mouseNo=1:length(LFP.ch(chNo).bandwidth(ii_plot).mouse)
                if group(mouseNo)==grNo
                    mouse_here=1;
                    this_mean_dB_power=zeros(1,length(LFP.ch(chNo).bandwidth(ii_plot).mouse(mouseNo).mean_dB_power));
                    this_mean_dB_power(1,:)=LFP.ch(chNo).bandwidth(ii_plot).mouse(mouseNo).mean_dB_power;

                    conv_this_mean_dB_power=conv(this_mean_dB_power,v,'same');

                    treatment_start=FileNames.mouse(mouseNo).treatment_start(1);
                    %Zero dB by subtracting the average of the day before
                    t_zero_start=LFP.mice(mouseNo).t_night_off(floor(treatment_start/24)-1);
                    t_zero_end=LFP.mice(mouseNo).t_night_on(floor(treatment_start/24));
                    dB_zero=mean(conv_this_mean_dB_power((t_hours>=t_zero_start)&(t_hours<=t_zero_end)));
                    if isnan(dB_zero)
                        dB_zero=0;
                    end
                    conv_this_mean_dB_power=conv_this_mean_dB_power-dB_zero;
                    skip_conv_this_mean_dB_power=conv_this_mean_dB_power(skip_ii+1:end-skip_ii);

                    if override_lims==0
                        min_dB=prctile(skip_conv_this_mean_dB_power,this_per)-0.1*(prctile(skip_conv_this_mean_dB_power,100-this_per)-prctile(skip_conv_this_mean_dB_power,this_per));
                        max_dB=prctile(skip_conv_this_mean_dB_power,100-this_per)+0.1*(prctile(skip_conv_this_mean_dB_power,100-this_per)-prctile(skip_conv_this_mean_dB_power,this_per));
                    else
                        min_dB=fig_ylim(1);
                        max_dB=fig_ylim(2);
                    end

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
                    for ts_ii=1:length(FileNames.mouse(mouseNo).treatment_start)
                        treatment_start=FileNames.mouse(mouseNo).treatment_start(ts_ii);
                        plot([treatment_start treatment_start], [min_dB max_dB],'-r','LineWidth',3)
                    end


                end

            end


         
            title(group_name{grNo})

            if mouse_here==1
                ylim([min_dB max_dB])
                if override_lims==1
                    xlim([fig_xlim])
                end
                ylabel('dB')
                xlabel('time (hrs)')
            end
        end


    end
end

%Plot delta dB bar graphs


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

        sgtitle(['Delta night/day power for ' einfSignalLabels{chNo} ' ' f_label{ii_plot}])

        for grNo=1:max(group)

            subplot(max(group),1,grNo)
            hold on

            mouse_here=0;
            for mouseNo=1:length(LFP.ch(chNo).bandwidth(ii_plot).mouse)
                if group(mouseNo)==grNo
                    
                    for ii_delta=1:length(LFP.ch(chNo).bandwidth(ii_plot).mouse(mouseNo).deltadBtimepoint)
                        bar(ii_delta+1,LFP.ch(chNo).bandwidth(ii_plot).mouse(mouseNo).deltadBtimepoint(ii_delta).delta_dB_power_night_day,'b')
                    end
                    
                    
                    %Place a vertical bar at start of treatment
                    for ts_ii=1:length(FileNames.mouse(mouseNo).treatment_start)
                        treatment_start=FileNames.mouse(mouseNo).treatment_start(ts_ii)/24 +0.2;
                        plot([treatment_start treatment_start], [-4 4],'-r','LineWidth',3)
                    end


                end

            end


         
            title(group_name{grNo})

            if mouse_here==1
                ylim([-4 4])
           
                ylabel('delta dB')
                xlabel('night number')
            end
        end


    end
end

%Plot normalized delta dB
%Plot delta dB bar graphs


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

        sgtitle(['Normalized delta night/day power for ' einfSignalLabels{chNo} ' ' f_label{ii_plot}])

        for grNo=1:max(group)

            subplot(max(group),1,grNo)
            hold on

            mouse_here=0;
            for mouseNo=1:length(LFP.ch(chNo).bandwidth(ii_plot).mouse)
                if group(mouseNo)==grNo
                    
                    all_deltas=[];
                    for ii_delta=1:length(LFP.ch(chNo).bandwidth(ii_plot).mouse(mouseNo).deltadBtimepoint)
                        all_deltas=[all_deltas LFP.ch(chNo).bandwidth(ii_plot).mouse(mouseNo).deltadBtimepoint(ii_delta).delta_dB_power_night_day];
                        
                    end
                    all_deltas=all_deltas/mean(all_deltas(3:5));
                    for ii_delta=1:length(LFP.ch(chNo).bandwidth(ii_plot).mouse(mouseNo).deltadBtimepoint)
                        bar(ii_delta+1,all_deltas(ii_delta),'b')
                    end

                    %Place a vertical bar at start of treatment
                    for ts_ii=1:length(FileNames.mouse(mouseNo).treatment_start)
                        treatment_start=FileNames.mouse(mouseNo).treatment_start(ts_ii)/24 +0.2;
                        plot([treatment_start treatment_start], [0 2],'-r','LineWidth',3)
                    end


                end

            end


         
            title(group_name{grNo})

            if mouse_here==1
                ylim([0 2])
           
                ylabel('Normalized delta dB')
                xlabel('night number')
            end
        end


    end
end

fprintf(1, 'The LFP data were convolved within a window of %d hours\n',conv_window*(t_hours(2)-t_hours(1)));

pffft=1;
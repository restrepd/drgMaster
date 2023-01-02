%drgMultiDayLFP

PathName='/Users/restrepd/Documents/Projects/Multi day LFP/';

figNo=0;

skip_ii=250; %Skip the first and last 10 sec

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
f_label{4}='Gamma_low';
f_label{5}='Gamma_high';

FileName{1}='Restrepo C031_2022-10-26_10_00_14_export_ov0p5_fft.mat';
group(1)=1;
treatment_start(1)=85;

FileName{2}='Restrepo T034_2022-10-26_10_00_29_export_ov0p5_fft.mat';
group(2)=2;
treatment_start(2)=85;

FileName{3}='Restrepo T035_2022-10-26_10_00_19_export_ov0p5_fft.mat';
group(3)=3;
treatment_start(3)=85;

LFPs=[];

%Plot dB per bandwidth
%The dB timecourse will be convolved with a flat window
for fileNo=1:length(FileName)

    load([PathName FileName{fileNo}])
    t_hours=handlespf.drgbchoices.start_hour+[1:size(mean_dB_power,2)]*(9/(60*60));
    LFP.files(fileNo).t_hours=t_hours;
    t=handlespf.drgbchoices.light_on;

    %This only works if the start is during the day
    LFP.files(fileNo).t_night_on=[];
    LFP.files(fileNo).t_night_off=[];

    ii_night=0;
    while t< t_hours(end)

        if (rem(t,24)<handlespf.drgbchoices.light_on)||(rem(t,24)>=handlespf.drgbchoices.light_off)

            %This is the night
            delta_t_night_end=24-handlespf.drgbchoices.light_off+handlespf.drgbchoices.light_on;
            ii_night=ii_night+1;
            if t+delta_t_night_end<t_hours(end)
                LFP.files(fileNo).t_night_on(ii_night)=t;
                LFP.files(fileNo).t_night_off(ii_night)=t+delta_t_night_end;
            else
                LFP.files(fileNo).t_night_on(ii_night)=t;
                LFP.files(fileNo).t_night_off(ii_night)=t_hours(end);
            end

            t=t+delta_t_night_end;

        else
            %This is the day
            t=t+(handlespf.drgbchoices.light_off-rem(t,24));
        end

    end

    for chNo=1:4

        for ii_plot=1:5
            this_mean_dB_power=zeros(1,size(mean_dB_power,2));
            this_mean_dB_power(1,:)=mean_dB_power(ii_plot,:,chNo);
            LFP.ch(chNo).bandwidth(ii_plot).file(fileNo).mean_dB_power=this_mean_dB_power;
        end
    end

end

this_per=5;
v=ones(1,500)/500;
for chNo=1:4

    for ii_plot=1:5

        this_mean_dB_power=[];
        for fileNo=1:length(LFP.ch(chNo).bandwidth(ii_plot).file)
            conv_dB=conv(LFP.ch(chNo).bandwidth(ii_plot).file(fileNo).mean_dB_power(skip_ii+1:end-skip_ii),v,'same');
            this_mean_dB_power=[this_mean_dB_power conv_dB];
        end

        LFP.ch(chNo).bandwidth(ii_plot).min_dB=prctile(this_mean_dB_power,this_per)-0.1*(prctile(this_mean_dB_power,100-this_per)-prctile(this_mean_dB_power,this_per));
        LFP.ch(chNo).bandwidth(ii_plot).max_dB=prctile(this_mean_dB_power,100-this_per)+0.1*(prctile(this_mean_dB_power,100-this_per)-prctile(this_mean_dB_power,this_per));
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
        set(hFig, 'units','normalized','position',[.07 .05 .45 .75])
        hold on

        sgtitle(['Power for ' einfSignalLabels{chNo} ' ' f_label{ii_plot}])

        for grNo=1:max(group)

            subplot(max(group),1,grNo)
            hold on
            v=ones(1,500)/500;

            for fileNo=1:length(LFP.ch(chNo).bandwidth(ii_plot).file)
                if group(fileNo)==grNo
                    this_mean_dB_power=zeros(1,length(LFP.ch(chNo).bandwidth(ii_plot).file(fileNo).mean_dB_power));
                    this_mean_dB_power(1,:)=LFP.ch(chNo).bandwidth(ii_plot).file(fileNo).mean_dB_power;

                    conv_this_mean_dB_power=conv(this_mean_dB_power,v,'same');

                    %Align the time to the start of the first night
                    t_hours=LFP.files(fileNo).t_hours-LFP.files(fileNo).t_night_on(1);
                    plot(t_hours(skip_ii+1:end-skip_ii),conv_this_mean_dB_power(skip_ii+1:end-skip_ii),'-k')
                end

            end

            min_dB=LFP.ch(chNo).bandwidth(ii_plot).min_dB;
            max_dB=LFP.ch(chNo).bandwidth(ii_plot).max_dB;

            for ii_nights=1:length(LFP.files(1).t_night_on)
                plot([LFP.files(1).t_night_on(ii_nights)-LFP.files(fileNo).t_night_on(1) LFP.files(1).t_night_off(ii_nights)-LFP.files(fileNo).t_night_on(1)],...
                    [min_dB+0.1*(prctile(this_mean_dB_power,99)-prctile(this_mean_dB_power,1))...
                    min_dB+0.1*(prctile(this_mean_dB_power,99)-prctile(this_mean_dB_power,1)) ], '-k', 'LineWidth', 3)
            end
            title(group_name{grNo})

            ylim([min_dB max_dB])
            ylabel('dB')
            xlabel('time (hrs)')
        end


    end
end
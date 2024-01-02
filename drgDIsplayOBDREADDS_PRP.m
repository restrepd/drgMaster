function drgDIsplayOBDREADDS_PRP
%drgOBChr2LFPpower 
%provides power analysis for Joe's ChR2 runs

clear all

[FileName,PathName] = uigetfile({'prp_*.mat'},'Select the drgOBDREADDS_PRP prp_ .mat file ');

load([PathName FileName])

ii_start_per_file=handles_out.ii_start_per_file;
ii_end_per_file=handles_out.ii_end_per_file;
time=handles_out.time;

close all
figNo=0;

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

plot(time(ii_start_per_file(1):ii_end_per_file(1)),mod_indx(ii_start_per_file(1):ii_end_per_file(1)),'ob-')
plot(time(ii_start_per_file(2):ii_end_per_file(2)),mod_indx(ii_start_per_file(2):ii_end_per_file(2)),'ob-')

if length(mod_indx)>20
    no_conv_points=6;
    conv_win=ones(1,no_conv_points);
    mod_indx_extend=[mean(mod_indx(1:no_conv_points/2))*ones(1,no_conv_points) mod_indx mean(mod_indx(end-(no_conv_points/2)+1:end))*ones(1,no_conv_points)];
    conv_mod_indx_extend=conv(mod_indx_extend,conv_win)/no_conv_points;
    conv_mod_indx=conv_mod_indx_extend(no_conv_points+1:no_conv_points+length(mod_indx));
    plot(time(ii_start_per_file(1):ii_end_per_file(1)),conv_mod_indx(ii_start_per_file(1):ii_end_per_file(1)),'-b','LineWidth',3)
    plot(time(ii_start_per_file(2):ii_end_per_file(2)),conv_mod_indx(ii_start_per_file(2):ii_end_per_file(2)),'-b','LineWidth',3)
end
theseyl=ylim;
between_ii=(time(ii_end_per_file(1))+time(ii_start_per_file(2)))/2;
plot([between_ii between_ii],theseyl,'-k','LineWidth',3)
xlabel('Time (sec)')
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

reference_power=(mean(handles_out.mean_troughPower(ii_start_per_file(1):ii_end_per_file(1)))+...
    mean(handles_out.mean_peakPower(ii_start_per_file(1):ii_end_per_file(1))))/2;
meanTroughPower=handles_out.mean_troughPower-reference_power;
meanPeakPower=handles_out.mean_peakPower-reference_power;

plot(time(ii_start_per_file(1):ii_end_per_file(1)),meanTroughPower(ii_start_per_file(1):ii_end_per_file(1)),'b-')
plot(time(ii_start_per_file(2):ii_end_per_file(2)),meanTroughPower(ii_start_per_file(2):ii_end_per_file(2)),'b-')

plot(time(ii_start_per_file(1):ii_end_per_file(1)),meanPeakPower(ii_start_per_file(1):ii_end_per_file(1)),'r-')
plot(time(ii_start_per_file(2):ii_end_per_file(2)),meanPeakPower(ii_start_per_file(2):ii_end_per_file(2)),'r-')

plot([between_ii between_ii],theseyl,'-k','LineWidth',3)

this_xl=xlim;
plot([this_xl],[0 0],'-k')

xlabel('Time(sec)')
ylabel('Peak power (dB)')


legend('Through','Peak')
title(['Power (dB), Do not use'])

%Peak and trough processed through PAC
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig2 = figure(figNo);
set(hFig2, 'units','normalized','position',[.25 .1 .4 .25])
hold on

reference_power=(mean(handles_out.troughPowerPAC(ii_start_per_file(1):ii_end_per_file(1)))+...
    mean(handles_out.peakPowerPAC(ii_start_per_file(1):ii_end_per_file(1))))/2;
meanTroughPower=handles_out.troughPowerPAC-reference_power;
meanPeakPower=handles_out.peakPowerPAC-reference_power;

plot(time(ii_start_per_file(1):ii_end_per_file(1)),meanTroughPower(ii_start_per_file(1):ii_end_per_file(1)),'b-')
plot(time(ii_start_per_file(2):ii_end_per_file(2)),meanTroughPower(ii_start_per_file(2):ii_end_per_file(2)),'b-')

plot(time(ii_start_per_file(1):ii_end_per_file(1)),meanPeakPower(ii_start_per_file(1):ii_end_per_file(1)),'r-')
plot(time(ii_start_per_file(2):ii_end_per_file(2)),meanPeakPower(ii_start_per_file(2):ii_end_per_file(2)),'r-')

plot([between_ii between_ii],theseyl,'-k','LineWidth',3)

this_xl=xlim;
plot([this_xl],[0 0],'-k')

xlabel('Time(sec)')
ylabel('Peak power (dB)')


legend('Through','Peak')
title(['Power (dB, PAC calculated)  '])


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
plot(time(ii_start_per_file(2):ii_end_per_file(2)),troughAngle(ii_start_per_file(2):ii_end_per_file(2)),'bo-')

plot(time(ii_start_per_file(1):ii_end_per_file(1)),troughAngle_conv(ii_start_per_file(1):ii_end_per_file(1)),'b-','LineWidth',3)
plot(time(ii_start_per_file(2):ii_end_per_file(2)),troughAngle_conv(ii_start_per_file(2):ii_end_per_file(2)),'b-','LineWidth',3)

plot(time(ii_start_per_file(1):ii_end_per_file(1)),peakAngle(ii_start_per_file(1):ii_end_per_file(1)),'ro-')
plot(time(ii_start_per_file(2):ii_end_per_file(2)),peakAngle(ii_start_per_file(2):ii_end_per_file(2)),'ro-')

plot(time(ii_start_per_file(1):ii_end_per_file(1)),peakAngle_conv(ii_start_per_file(1):ii_end_per_file(1)),'r-','LineWidth',3)
plot(time(ii_start_per_file(2):ii_end_per_file(2)),peakAngle_conv(ii_start_per_file(2):ii_end_per_file(2)),'r-','LineWidth',3)

plot([between_ii between_ii],theseyl,'-k','LineWidth',3)

xlabel('Time(sec)')
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
xlabel('Time (sec)')
ylabel('Phase for low freq oscillation (deg)');
title(['Phase-amplitude coupling timecourse '])

pffft=1
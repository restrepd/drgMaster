function drgLickTimecourse(handles)
%Performs an event-related analysis. The event is signaled by a sharp chane
%in the reference voltage. This is used to analyze lick-related changes in
%LFP

% 
% evTypeNo=handles.evTypeNo;
% [per_corr_per_trialEv1, lick_timecourseEv1,which_eventEv1,no_trialsEv1,inter_lick_intervals]=drgLickTimecourseThisEv(handles);
% 
% handles.evTypeNo=handles.evTypeNo2;
% [per_corr_per_trialEv2, lick_timecourseEv2,which_eventEv2,no_trialsEv2,inter_lick_intervals]=drgLickTimecourseThisEv(handles);
% 
% handles.evTypeNo=evTypeNo;
% 
% %Timecourse for licks
% try
%     close 4
% catch
% end
% 
% hFig4 = figure(4);
% set(hFig4, 'units','normalized','position',[.05 .4 .5 .25])
% 
% no_lick_dt=floor(((handles.time_end-handles.time_pad)-(handles.time_start+handles.time_pad))/handles.dt_lick);
% time=[1:no_lick_dt]*handles.dt_lick+handles.time_start+handles.time_pad;
% 
% Ev1_mean_licks=mean(lick_timecourseEv1,1)/handles.dt_lick;
% Ev1_CI_licks = bootci(1000, @mean, lick_timecourseEv1)/handles.dt_lick;
% mean_Ev1=[Ev1_mean_licks;Ev1_mean_licks];
% Ev1_ci_l=Ev1_CI_licks(1,:)-mean_Ev1;
% 
% [hl1, hp1] = boundedline(time',Ev1_mean_licks', Ev1_ci_l', 'r');
% 
% hold on
% 
% Ev2_mean_licks=mean(lick_timecourseEv2,1)/handles.dt_lick;
% Ev2_CI_licks = bootci(1000, @mean, lick_timecourseEv2)/handles.dt_lick;
% mean_Ev2=[Ev2_mean_licks;Ev2_mean_licks];
% Ev2_ci_l=Ev2_CI_licks(1,:)-mean_Ev2;
% 
% [hl2, hp2] = boundedline(time',Ev2_mean_licks', Ev2_ci_l', 'b');
% 
% title('Lick timecourse, S+ red, S- blue')
% xlabel('Time (s)')
% ylabel('Licks/s')
% set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)

evTypeNo=handles.evTypeNo;
[lick_freq1,times_lick_freq,lick_traces1,CIlickf1,lick_trace_times,stamped_lick_ii,these_stamped_lick_times,no_trials,trials_included]=drgGetLicks(handles);


handles.evTypeNo=handles.evTypeNo2;
[lick_freq2,times_lick_freq,lick_traces2,CIlickf2,lick_trace_times,stamped_lick_ii,these_stamped_lick_times,no_trials,trials_included]=drgGetLicks(handles);

handles.evTypeNo=evTypeNo;

%Timecourse for licks
try
    close 4
catch
end

hFig4 = figure(4);
set(hFig4, 'units','normalized','position',[.05 .4 .5 .25])


%plot(times_lick_freq, lick_freq)
[hl1, hp1] = boundedline(times_lick_freq',lick_freq1', CIlickf1', 'r');

[hl1, hp1] = boundedline(times_lick_freq',lick_freq2', CIlickf2', 'b');

ylim([0 1.2*max([max(lick_freq1) max(lick_freq1)])])
xlabel('Time (sec)')
ylabel('frequency (Hz)')
title('Lick frequency')

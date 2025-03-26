function drgOBtoHip_summary_batch_exclude_tPRP_V3(choiceBatchPathName,choiceFileName)

%Note: fitcnet will not work in Matlab versions earlier than 2021a

close all

clear all

PathName_out='/Users/restrepd/Documents/Projects/JoeV/2025 3 VZV 2 Mock/';

first_file=1;

LFPPowerFreqLow=[6 15 35 65];

LFPPowerFreqHigh=[14 30 55 95];

badwidth_names{1}='Theta';

badwidth_names{2}='Beta';

badwidth_names{3}='Low gamma';

badwidth_names{4}='High gamma';

actual_electrode(1)=1;

actual_electrode(2)=8;

actual_electrode(3)=9;

actual_electrode(4)=16;

our_colors(1).color=[0.9 0.6 0];

our_colors(2).color=[0.35 0.7 0.9];

our_colors(3).color=[0 0.6 0.5];

our_colors(4).color=[0.95 0.9 0.25];

%Re-arrange the groups

%Group numbers

%1 WT Treated

%2 5xFAD Treated

%3 WT Untreated

%4 5xFAD Untreated

pre_post_label{1}='Pre';

pre_post_label{2}='Post';

week_group_label={'Exo Weeks 1&2','Exo Week 3','Poly dA:dT'}

group_label{1}='Mock';

group_label{2}='VZV';

group_label{3}='';

group_label{4}='';

genotype_label='WT';

%genotype_label{2}='WT';

laser_label='T';

%laser_label{2}='T';

time_window_labels = {'Pre','Laser', 'Post'}; % Add all your group labels here

logP_time_window_labels = {'Pre','Laser', 'Post'};

t_ref_start=0;

t_ref_end=20*60; %In seconds

t_laser_end=40*60;

%You can choose which plots are shown

these_plots_shown(1)=0; %delta peak tPRP bar graph for laser and post windows (minus pre)

these_plots_shown(2)=0; %delta peak tPRP bar graph for laser window (minus pre) per week

these_plots_shown(3)=0; %delta trough tPRP bar graph for laser and post windows (minus pre)

these_plots_shown(4)=0; %delta trough tPRP bar graph for laser window (minus pre)

these_plots_shown(5)=1; %tPRP peak for each electrode and pre- and laser windows per week

these_plots_shown(6)=0; %tPRP peak for each electrode for laser window only

these_plots_shown(7)=1; %tPRP trough for each electrode and pre- and laser window per week

these_plots_shown(8)=0; %tPRP trough for each electrode for laser window only

these_plots_shown(9)=1; %show modulation index for each electrode for each time window

these_plots_shown(10)=0; %show modulation index for each electrode for each time window

these_plots_shown(11)=0; %show modulation index as a function of session number for each electrode for only one time window

these_plots_shown(12)=1; %show peak angle variance for each electrode for each time window

these_plots_shown(13)=0; %show peak angle variance as a function of session number for each electrode for only one time window

these_plots_shown(14)=1; %tPRP peak for pre-week 1 and post-week 4 only line 1602

these_plots_shown(15)=1; %tPRP trough for pre-week 1 and post-week 4 only

% Define trials (weeks)
trials = {'pre-week1', 'pre-week2', 'pre-week3', 'pre-week4', 'post-week1', 'post-week2', 'post-week3', 'post-week4'};


%Diego has the choices file in /Users/restrepd/Documents/Projects/Joe_OB_to_hippo/5xFADvsWT_1_hour_treatmentVsnone

if nargin==0

[choiceFileName,choiceBatchPathName] = uigetfile({'drgbChoicesOBtoHIP*.m'},'Select the .m file with all the choices for analysis');

end

fprintf(1, ['\ndrgOBtoHIP_batch_tPRP run for ' choiceFileName '\n\n']);

addpath(choiceBatchPathName)

eval(['handles=' choiceFileName(1:end-2) ';'])

handles.choiceFileName=choiceFileName;

handles.choiceBatchPathName=choiceBatchPathName;

new_no_files=handles.no_files;

%Parallel batch processing for each file

all_files_present=1;

for filNum=first_file:handles.no_files

%Make sure that all the files exist

jtFileName=handles.FileName{filNum};

this_filename=handles.FileName{filNum};

FileName_out=['OBtoHIPtPRP_' this_filename(11:end)]

if exist([PathName_out FileName_out])==0

fprintf(1, ['Program will be terminated because file No %d, does not exist\n'],filNum);

all_files_present=0;

end

end

%Text file for statistical output

fileID = fopen([choiceBatchPathName 'OBtoHIP_summary_batch_LFP_stats.txt'],'w');

tic

if all_files_present==1

these_mice=unique(handles.mouse_no);

figNo=0;

groups=unique(handles.group_no);


%--------------------------------------------------------------------------
% Graph 1 & 2: Delta Peak/Trough tPRP for HIP
%--------------------------------------------------------------------------

%Do hippocampus
ii_electrodes=[1 2]; %This HIP

%Peak
figNo = figNo + 1;
hFig = figure(figNo);
set(hFig, 'units', 'normalized', 'position', [.1 .1 .4 .4]); % Adjust figure size
hold on

title(['Delta Peak tPRP by Trial - Hippocampus']);
ylabel('Delta tPRP (dB)');

these_trials(1).trials=[1 2];
these_trials(2).trials=[5 6];
these_trials(3).trials=3;
these_trials(4).trials=7;
these_trials(5).trials=4;
these_trials(6).trials=8;

x_pos=0;

glm_tPRP=[];
glm_tPRP_ii=0;

tPRP_data=[];
ii_tPRP=0;

for ii=1:2:5
    week_group=floor(ii/2)+1;
    for groupNo = 1:length(groups)
        % Collect data for this group and trial
        these_delta_tPRP_pre = [];
        these_delta_tPRP_post = [];

        for ii_indx=ii:ii+1
            for fileNo = 1:handles.no_files
                if (handles.group_no(fileNo) == groupNo)&(sum(handles.session_no(fileNo) == these_trials(ii_indx).trials)>0)

                    %  %Check if this file contains the appropriate week number
                    %  file_week_number = handles.file(fileNo).week_number; %Assumes you have a week_number field in your handles.file struct
                    %
                    % %Skip files from the wrong week
                    %  if file_week_number ~= week_number
                    %      continue;
                    %  end

                    this_filename = handles.FileName{fileNo};
                    FileName_out = ['OBtoHIPtPRP_' this_filename(11:end)];
 
                    % if exist([PathName_out FileName_out], 'file') %Take out
                    load([PathName_out FileName_out]);

                    for ii_electrode=ii_electrodes
                        add_electrode = (sum(handles.file(fileNo).exclude_electrodes==actual_electrode(ii_electrode))==0);

                        if add_electrode
                            log_P_timecourse = handles_out.electrode(ii_electrode).handles_out.mean_peakPower;
                            time = handles_out.electrode(ii_electrode).handles_out.time;

                            % Calculate delta tPRP (Laser - Pre)
                            delta_tPRP = mean(log_P_timecourse((time >= t_ref_end) & (time <= t_laser_end))) - ...
                                mean(log_P_timecourse(time <= t_ref_end));

                            if (sum(handles.session_no(fileNo) == these_trials(ii).trials)>0)
                                %Pre
                                these_delta_tPRP_pre = [these_delta_tPRP_pre, delta_tPRP];
                            else
                                %Post
                                these_delta_tPRP_post = [these_delta_tPRP_post, delta_tPRP];
                            end
                        end
                    end
                    % end
                end
            end
        end

        % Plotting
        bar(x_pos, mean(these_delta_tPRP_pre), 0.8, 'FaceColor', our_colors(groupNo).color, 'EdgeColor', 'none'); % Use group color
        CI_points=bootci(1000, @mean, these_delta_tPRP_pre);
        plot([x_pos x_pos],CI_points,'-k','LineWidth',3)

        x_pos=x_pos+1;


        bar(x_pos, mean(these_delta_tPRP_post), 0.8, 'FaceColor', our_colors(groupNo).color, 'EdgeColor', 'none'); % Use group color
        CI_points=bootci(1000, @mean, these_delta_tPRP_post);
        plot([x_pos x_pos],CI_points,'-k','LineWidth',3)


        for jj=1:length(these_delta_tPRP_pre)
            plot(x_pos-1,these_delta_tPRP_pre(jj),'ok','MarkerSize',5,'MarkerFaceColor','k')
            plot(x_pos,these_delta_tPRP_post(jj),'ok','MarkerSize',5,'MarkerFaceColor','k')
            plot([x_pos-1 x_pos],[these_delta_tPRP_pre(jj) these_delta_tPRP_post(jj)],'-k','LineWidth',1)
        end

        x_pos=x_pos+1;

        %Enter data for glm
        glm_tPRP.data(glm_tPRP_ii+1:glm_tPRP_ii+length(these_delta_tPRP_pre))=these_delta_tPRP_pre;
        glm_tPRP.group(glm_tPRP_ii+1:glm_tPRP_ii+length(these_delta_tPRP_pre))=groupNo;
        glm_tPRP.pre_post(glm_tPRP_ii+1:glm_tPRP_ii+length(these_delta_tPRP_pre))=1;
        glm_tPRP.week(glm_tPRP_ii+1:glm_tPRP_ii+length(these_delta_tPRP_pre))=week_group;
        glm_tPRP_ii=glm_tPRP_ii+length(these_delta_tPRP_pre);

        glm_tPRP.data(glm_tPRP_ii+1:glm_tPRP_ii+length(these_delta_tPRP_post))=these_delta_tPRP_post;
        glm_tPRP.group(glm_tPRP_ii+1:glm_tPRP_ii+length(these_delta_tPRP_post))=groupNo;
        glm_tPRP.pre_post(glm_tPRP_ii+1:glm_tPRP_ii+length(these_delta_tPRP_post))=2;
        glm_tPRP.week(glm_tPRP_ii+1:glm_tPRP_ii+length(these_delta_tPRP_post))=week_group;
        glm_tPRP_ii=glm_tPRP_ii+length(these_delta_tPRP_post);

        %Enter the data for t-test/ranksum
        ii_tPRP=ii_tPRP+1;
        tPRP_data(ii_tPRP).data=these_delta_tPRP_pre;
        tPRP_data(ii_tPRP).description=[group_label{groupNo} ' Pre  ' week_group_label{week_group}];

        ii_tPRP=ii_tPRP+1;
        tPRP_data(ii_tPRP).data=these_delta_tPRP_post;
        tPRP_data(ii_tPRP).description=[group_label{groupNo} ' Post  ' week_group_label{week_group}];            

    end
    x_pos=x_pos+1;
end
xticks([0 1 2 3 5 6 7 8 10 11 12 13]);
xticklabels({'Pre', 'Post', 'Pre', 'Post', 'Pre', 'Post', 'Pre', 'Post',...
    'Pre', 'Post', 'Pre', 'Post' });
ylim([-8 7])
text(-0.5, 5.5, 'Exo Weeks 1&2', 'FontSize', 14, 'FontWeight', 'bold')
text(5, 5.5, 'Exo Week 3', 'FontSize', 14, 'FontWeight', 'bold')
text(9.5, 5.5, 'Poly dA:dT Week 4', 'FontSize', 14, 'FontWeight', 'bold')
text(11,-4,'Mock', 'FontSize', 14, 'FontWeight', 'bold', 'Color',our_colors(1).color)
text(11,-5,'VZV', 'FontSize', 14, 'FontWeight', 'bold', 'Color',our_colors(2).color)

%Now do the glm
fprintf(1, ['\n\nglm for peak delta tPRP hippocampus\n'])
tbl = table(glm_tPRP.data',glm_tPRP.group',glm_tPRP.pre_post',glm_tPRP.week',...
    'VariableNames',{'tPRP','group','pre_post','week'});
mdl = fitglm(tbl,'tPRP~group+pre_post+week+group*pre_post*week'...
    ,'CategoricalVars',[2,3,4])


fprintf(1, ['\n\nRanksum or t-test for peak delta tPRP hippocampus \n'])
%Now do the ranksums
output_data = drgMutiRanksumorTtest(tPRP_data);

%Trough
figNo = figNo + 1;
hFig = figure(figNo);
set(hFig, 'units', 'normalized', 'position', [.1 .1 .4 .4]); % Adjust figure size
hold on

title(['Delta Trough tPRP by Trial - Hippocampus']);
ylabel('Delta tPRP (dB)');

these_trials(1).trials=[1 2];
these_trials(2).trials=[5 6];
these_trials(3).trials=3;
these_trials(4).trials=7;
these_trials(5).trials=4;
these_trials(6).trials=8;

x_pos=0;

glm_tPRP=[];
glm_tPRP_ii=0;

tPRP_data=[];
ii_tPRP=0;

for ii=1:2:5
    week_group=floor(ii/2)+1;
    for groupNo = 1:length(groups)
        % Collect data for this group and trial
        these_delta_tPRP_pre = [];
        these_delta_tPRP_post = [];

        for ii_indx=ii:ii+1
            for fileNo = 1:handles.no_files
                if (handles.group_no(fileNo) == groupNo)&(sum(handles.session_no(fileNo) == these_trials(ii_indx).trials)>0)

                    %  %Check if this file contains the appropriate week number
                    %  file_week_number = handles.file(fileNo).week_number; %Assumes you have a week_number field in your handles.file struct
                    %
                    % %Skip files from the wrong week
                    %  if file_week_number ~= week_number
                    %      continue;
                    %  end

                    this_filename = handles.FileName{fileNo};
                    FileName_out = ['OBtoHIPtPRP_' this_filename(11:end)];
 
                    % if exist([PathName_out FileName_out], 'file') %Take out
                    load([PathName_out FileName_out]);

                    for ii_electrode=ii_electrodes
                        add_electrode = (sum(handles.file(fileNo).exclude_electrodes==actual_electrode(ii_electrode))==0);

                        if add_electrode
                            log_P_timecourse = handles_out.electrode(ii_electrode).handles_out.mean_troughPower;
                            time = handles_out.electrode(ii_electrode).handles_out.time;

                            % Calculate delta tPRP (Laser - Pre)
                            delta_tPRP = mean(log_P_timecourse((time >= t_ref_end) & (time <= t_laser_end))) - ...
                                mean(log_P_timecourse(time <= t_ref_end));

                            if (sum(handles.session_no(fileNo) == these_trials(ii).trials)>0)
                                %Pre
                                these_delta_tPRP_pre = [these_delta_tPRP_pre, delta_tPRP];
                            else
                                %Post
                                these_delta_tPRP_post = [these_delta_tPRP_post, delta_tPRP];
                            end
                        end
                    end
                    % end
                end
            end
        end

        % Plotting


        bar(x_pos, mean(these_delta_tPRP_pre), 0.8, 'FaceColor', our_colors(groupNo).color, 'EdgeColor', 'none'); % Use group color
        CI_points=bootci(1000, @mean, these_delta_tPRP_pre);
        plot([x_pos x_pos],CI_points,'-k','LineWidth',3)

        x_pos=x_pos+1;


        bar(x_pos, mean(these_delta_tPRP_post), 0.8, 'FaceColor', our_colors(groupNo).color, 'EdgeColor', 'none'); % Use group color
        CI_points=bootci(1000, @mean, these_delta_tPRP_post);
        plot([x_pos x_pos],CI_points,'-k','LineWidth',3)


        for jj=1:length(these_delta_tPRP_pre)
            plot(x_pos-1,these_delta_tPRP_pre(jj),'ok','MarkerSize',5,'MarkerFaceColor','k')
            plot(x_pos,these_delta_tPRP_post(jj),'ok','MarkerSize',5,'MarkerFaceColor','k')
            plot([x_pos-1 x_pos],[these_delta_tPRP_pre(jj) these_delta_tPRP_post(jj)],'-k','LineWidth',1)
        end

        x_pos=x_pos+1;

            %Enter data for glm
        glm_tPRP.data(glm_tPRP_ii+1:glm_tPRP_ii+length(these_delta_tPRP_pre))=these_delta_tPRP_pre;
        glm_tPRP.group(glm_tPRP_ii+1:glm_tPRP_ii+length(these_delta_tPRP_pre))=groupNo;
        glm_tPRP.pre_post(glm_tPRP_ii+1:glm_tPRP_ii+length(these_delta_tPRP_pre))=1;
        glm_tPRP.week(glm_tPRP_ii+1:glm_tPRP_ii+length(these_delta_tPRP_pre))=week_group;
        glm_tPRP_ii=glm_tPRP_ii+length(these_delta_tPRP_pre);

        glm_tPRP.data(glm_tPRP_ii+1:glm_tPRP_ii+length(these_delta_tPRP_post))=these_delta_tPRP_post;
        glm_tPRP.group(glm_tPRP_ii+1:glm_tPRP_ii+length(these_delta_tPRP_post))=groupNo;
        glm_tPRP.pre_post(glm_tPRP_ii+1:glm_tPRP_ii+length(these_delta_tPRP_post))=2;
        glm_tPRP.week(glm_tPRP_ii+1:glm_tPRP_ii+length(these_delta_tPRP_post))=week_group;
        glm_tPRP_ii=glm_tPRP_ii+length(these_delta_tPRP_post);

        %Enter the data for t-test/ranksum
        ii_tPRP=ii_tPRP+1;
        tPRP_data(ii_tPRP).data=these_delta_tPRP_pre;
        tPRP_data(ii_tPRP).description=[group_label{groupNo} ' Pre  ' week_group_label{week_group}];

        ii_tPRP=ii_tPRP+1;
        tPRP_data(ii_tPRP).data=these_delta_tPRP_post;
        tPRP_data(ii_tPRP).description=[group_label{groupNo} ' Post  ' week_group_label{week_group}]; 
    end
    x_pos=x_pos+1;
end

xticks([0 1 2 3 5 6 7 8 10 11 12 13]);
xticklabels({'Pre', 'Post', 'Pre', 'Post', 'Pre', 'Post', 'Pre', 'Post',...
    'Pre', 'Post', 'Pre', 'Post' });
ylim([-8 7])
text(-0.5, 5.5, 'Exo Weeks 1&2', 'FontSize', 14, 'FontWeight', 'bold')
text(5, 5.5, 'Exo Week 3', 'FontSize', 14, 'FontWeight', 'bold')
text(9.5, 5.5, 'Poly dA:dT Week 4', 'FontSize', 14, 'FontWeight', 'bold')
text(11,-4,'Mock', 'FontSize', 14, 'FontWeight', 'bold', 'Color',our_colors(1).color)
text(11,-5,'VZV', 'FontSize', 14, 'FontWeight', 'bold', 'Color',our_colors(2).color)

%Now do the glm
fprintf(1, ['\n\nglm for trough delta tPRP hippocampus\n'])
tbl = table(glm_tPRP.data',glm_tPRP.group',glm_tPRP.pre_post',glm_tPRP.week',...
    'VariableNames',{'tPRP','group','pre_post','week'});
mdl = fitglm(tbl,'tPRP~group+pre_post+week+group*pre_post*week'...
    ,'CategoricalVars',[2,3,4])


fprintf(1, ['\n\nRanksum or t-test for trough delta tPRP hippocampus \n'])
%Now do the ranksums
output_data = drgMutiRanksumorTtest(tPRP_data);

%Now do olfactory bulb
ii_electrodes=[3 4]; %This OB

%Peak
figNo = figNo + 1;
hFig = figure(figNo);
set(hFig, 'units', 'normalized', 'position', [.1 .1 .4 .4]); % Adjust figure size
hold on

title(['Delta Peak tPRP by Trial - OB']);
ylabel('Delta tPRP (dB)');

these_trials(1).trials=[1 2];
these_trials(2).trials=[5 6];
these_trials(3).trials=3;
these_trials(4).trials=7;
these_trials(5).trials=4;
these_trials(6).trials=8;

x_pos=0;

glm_tPRP=[];
glm_tPRP_ii=0;

tPRP_data=[];
ii_tPRP=0;

for ii=1:2:5
    week_group=floor(ii/2)+1;
    for groupNo = 1:length(groups)
        % Collect data for this group and trial
        these_delta_tPRP_pre = [];
        these_delta_tPRP_post = [];

        for ii_indx=ii:ii+1
            for fileNo = 1:handles.no_files
                if (handles.group_no(fileNo) == groupNo)&(sum(handles.session_no(fileNo) == these_trials(ii_indx).trials)>0)

                    %  %Check if this file contains the appropriate week number
                    %  file_week_number = handles.file(fileNo).week_number; %Assumes you have a week_number field in your handles.file struct
                    %
                    % %Skip files from the wrong week
                    %  if file_week_number ~= week_number
                    %      continue;
                    %  end

                    this_filename = handles.FileName{fileNo};
                    FileName_out = ['OBtoHIPtPRP_' this_filename(11:end)];
 
                    % if exist([PathName_out FileName_out], 'file') %Take out
                    load([PathName_out FileName_out]);

                    for ii_electrode=ii_electrodes
                        add_electrode = (sum(handles.file(fileNo).exclude_electrodes==actual_electrode(ii_electrode))==0);

                        if add_electrode
                            log_P_timecourse = handles_out.electrode(ii_electrode).handles_out.mean_peakPower;
                            time = handles_out.electrode(ii_electrode).handles_out.time;

                            % Calculate delta tPRP (Laser - Pre)
                            delta_tPRP = mean(log_P_timecourse((time >= t_ref_end) & (time <= t_laser_end))) - ...
                                mean(log_P_timecourse(time <= t_ref_end));

                            if (sum(handles.session_no(fileNo) == these_trials(ii).trials)>0)
                                %Pre
                                these_delta_tPRP_pre = [these_delta_tPRP_pre, delta_tPRP];
                            else
                                %Post
                                these_delta_tPRP_post = [these_delta_tPRP_post, delta_tPRP];
                            end
                        end
                    end
                    % end
                end
            end
        end

        % Plotting


        bar(x_pos, mean(these_delta_tPRP_pre), 0.8, 'FaceColor', our_colors(groupNo).color, 'EdgeColor', 'none'); % Use group color
        CI_points=bootci(1000, @mean, these_delta_tPRP_pre);
        plot([x_pos x_pos],CI_points,'-k','LineWidth',3)

        x_pos=x_pos+1;


        bar(x_pos, mean(these_delta_tPRP_post), 0.8, 'FaceColor', our_colors(groupNo).color, 'EdgeColor', 'none'); % Use group color
        CI_points=bootci(1000, @mean, these_delta_tPRP_post);
        plot([x_pos x_pos],CI_points,'-k','LineWidth',3)


        for jj=1:length(these_delta_tPRP_pre)
            plot(x_pos-1,these_delta_tPRP_pre(jj),'ok','MarkerSize',5,'MarkerFaceColor','k')
            plot(x_pos,these_delta_tPRP_post(jj),'ok','MarkerSize',5,'MarkerFaceColor','k')
            plot([x_pos-1 x_pos],[these_delta_tPRP_pre(jj) these_delta_tPRP_post(jj)],'-k','LineWidth',1)
        end

        x_pos=x_pos+1;

            %Enter data for glm
        glm_tPRP.data(glm_tPRP_ii+1:glm_tPRP_ii+length(these_delta_tPRP_pre))=these_delta_tPRP_pre;
        glm_tPRP.group(glm_tPRP_ii+1:glm_tPRP_ii+length(these_delta_tPRP_pre))=groupNo;
        glm_tPRP.pre_post(glm_tPRP_ii+1:glm_tPRP_ii+length(these_delta_tPRP_pre))=1;
        glm_tPRP.week(glm_tPRP_ii+1:glm_tPRP_ii+length(these_delta_tPRP_pre))=week_group;
        glm_tPRP_ii=glm_tPRP_ii+length(these_delta_tPRP_pre);

        glm_tPRP.data(glm_tPRP_ii+1:glm_tPRP_ii+length(these_delta_tPRP_post))=these_delta_tPRP_post;
        glm_tPRP.group(glm_tPRP_ii+1:glm_tPRP_ii+length(these_delta_tPRP_post))=groupNo;
        glm_tPRP.pre_post(glm_tPRP_ii+1:glm_tPRP_ii+length(these_delta_tPRP_post))=2;
        glm_tPRP.week(glm_tPRP_ii+1:glm_tPRP_ii+length(these_delta_tPRP_post))=week_group;
        glm_tPRP_ii=glm_tPRP_ii+length(these_delta_tPRP_post);

        %Enter the data for t-test/ranksum
        ii_tPRP=ii_tPRP+1;
        tPRP_data(ii_tPRP).data=these_delta_tPRP_pre;
        tPRP_data(ii_tPRP).description=[group_label{groupNo} ' Pre  ' week_group_label{week_group}];

        ii_tPRP=ii_tPRP+1;
        tPRP_data(ii_tPRP).data=these_delta_tPRP_post;
        tPRP_data(ii_tPRP).description=[group_label{groupNo} ' Post  ' week_group_label{week_group}]; 
    end
    x_pos=x_pos+1;
end
xticks([0 1 2 3 5 6 7 8 10 11 12 13]);
xticklabels({'Pre', 'Post', 'Pre', 'Post', 'Pre', 'Post', 'Pre', 'Post',...
    'Pre', 'Post', 'Pre', 'Post' });
ylim([-6 12])
text(-0.5, 10, 'Exo Weeks 1&2', 'FontSize', 14, 'FontWeight', 'bold')
text(5, 10, 'Exo Week 3', 'FontSize', 14, 'FontWeight', 'bold')
text(9.5, 10, 'Poly dA:dT Week 4', 'FontSize', 14, 'FontWeight', 'bold')
text(11,-2,'Mock', 'FontSize', 14, 'FontWeight', 'bold', 'Color',our_colors(1).color)
text(11,-3,'VZV', 'FontSize', 14, 'FontWeight', 'bold', 'Color',our_colors(2).color)

%Now do the glm
fprintf(1, ['\n\nglm for peak delta tPRP olfactory bulb\n'])
tbl = table(glm_tPRP.data',glm_tPRP.group',glm_tPRP.pre_post',glm_tPRP.week',...
    'VariableNames',{'tPRP','group','pre_post','week'});
mdl = fitglm(tbl,'tPRP~group+pre_post+week+group*pre_post*week'...
    ,'CategoricalVars',[2,3,4])


fprintf(1, ['\n\nRanksum or t-test for peak delta tPRP olfactory bulb \n'])

%Now do the ranksums
output_data = drgMutiRanksumorTtest(tPRP_data);

%Trough
figNo = figNo + 1;
hFig = figure(figNo);
set(hFig, 'units', 'normalized', 'position', [.1 .1 .4 .4]); % Adjust figure size
hold on

title(['Delta Trough tPRP by Trial - OB']);
ylabel('Delta tPRP (dB)');

these_trials(1).trials=[1 2];
these_trials(2).trials=[5 6];
these_trials(3).trials=3;
these_trials(4).trials=7;
these_trials(5).trials=4;
these_trials(6).trials=8;

x_pos=0;

glm_tPRP=[];
glm_tPRP_ii=0;

tPRP_data=[];
ii_tPRP=0;


for ii=1:2:5
    week_group=floor(ii/2)+1;
    for groupNo = 1:length(groups)
        % Collect data for this group and trial
        these_delta_tPRP_pre = [];
        these_delta_tPRP_post = [];

        for ii_indx=ii:ii+1
            for fileNo = 1:handles.no_files
                if (handles.group_no(fileNo) == groupNo)&(sum(handles.session_no(fileNo) == these_trials(ii_indx).trials)>0)

                    %  %Check if this file contains the appropriate week number
                    %  file_week_number = handles.file(fileNo).week_number; %Assumes you have a week_number field in your handles.file struct
                    %
                    % %Skip files from the wrong week
                    %  if file_week_number ~= week_number
                    %      continue;
                    %  end

                    this_filename = handles.FileName{fileNo};
                    FileName_out = ['OBtoHIPtPRP_' this_filename(11:end)];
 
                    % if exist([PathName_out FileName_out], 'file') %Take out
                    load([PathName_out FileName_out]);

                    for ii_electrode=ii_electrodes
                        add_electrode = (sum(handles.file(fileNo).exclude_electrodes==actual_electrode(ii_electrode))==0);

                        if add_electrode
                            log_P_timecourse = handles_out.electrode(ii_electrode).handles_out.mean_troughPower;
                            time = handles_out.electrode(ii_electrode).handles_out.time;

                            % Calculate delta tPRP (Laser - Pre)
                            delta_tPRP = mean(log_P_timecourse((time >= t_ref_end) & (time <= t_laser_end))) - ...
                                mean(log_P_timecourse(time <= t_ref_end));

                            if (sum(handles.session_no(fileNo) == these_trials(ii).trials)>0)
                                %Pre
                                these_delta_tPRP_pre = [these_delta_tPRP_pre, delta_tPRP];
                            else
                                %Post
                                these_delta_tPRP_post = [these_delta_tPRP_post, delta_tPRP];
                            end
                        end
                    end
                    % end
                end
            end
        end

        % Plotting


        bar(x_pos, mean(these_delta_tPRP_pre), 0.8, 'FaceColor', our_colors(groupNo).color, 'EdgeColor', 'none'); % Use group color
        CI_points=bootci(1000, @mean, these_delta_tPRP_pre);
        plot([x_pos x_pos],CI_points,'-k','LineWidth',3)

        x_pos=x_pos+1;


        bar(x_pos, mean(these_delta_tPRP_post), 0.8, 'FaceColor', our_colors(groupNo).color, 'EdgeColor', 'none'); % Use group color
        CI_points=bootci(1000, @mean, these_delta_tPRP_post);
        plot([x_pos x_pos],CI_points,'-k','LineWidth',3)


        for jj=1:length(these_delta_tPRP_pre)
            plot(x_pos-1,these_delta_tPRP_pre(jj),'ok','MarkerSize',5,'MarkerFaceColor','k')
            plot(x_pos,these_delta_tPRP_post(jj),'ok','MarkerSize',5,'MarkerFaceColor','k')
            plot([x_pos-1 x_pos],[these_delta_tPRP_pre(jj) these_delta_tPRP_post(jj)],'-k','LineWidth',1)
        end

        x_pos=x_pos+1;

        
            %Enter data for glm
        glm_tPRP.data(glm_tPRP_ii+1:glm_tPRP_ii+length(these_delta_tPRP_pre))=these_delta_tPRP_pre;
        glm_tPRP.group(glm_tPRP_ii+1:glm_tPRP_ii+length(these_delta_tPRP_pre))=groupNo;
        glm_tPRP.pre_post(glm_tPRP_ii+1:glm_tPRP_ii+length(these_delta_tPRP_pre))=1;
        glm_tPRP.week(glm_tPRP_ii+1:glm_tPRP_ii+length(these_delta_tPRP_pre))=week_group;
        glm_tPRP_ii=glm_tPRP_ii+length(these_delta_tPRP_pre);

        glm_tPRP.data(glm_tPRP_ii+1:glm_tPRP_ii+length(these_delta_tPRP_post))=these_delta_tPRP_post;
        glm_tPRP.group(glm_tPRP_ii+1:glm_tPRP_ii+length(these_delta_tPRP_post))=groupNo;
        glm_tPRP.pre_post(glm_tPRP_ii+1:glm_tPRP_ii+length(these_delta_tPRP_post))=2;
        glm_tPRP.week(glm_tPRP_ii+1:glm_tPRP_ii+length(these_delta_tPRP_post))=week_group;
        glm_tPRP_ii=glm_tPRP_ii+length(these_delta_tPRP_post);

        %Enter the data for t-test/ranksum
        ii_tPRP=ii_tPRP+1;
        tPRP_data(ii_tPRP).data=these_delta_tPRP_pre;
        tPRP_data(ii_tPRP).description=[group_label{groupNo} ' Pre  ' week_group_label{week_group}];

        ii_tPRP=ii_tPRP+1;
        tPRP_data(ii_tPRP).data=these_delta_tPRP_post;
        tPRP_data(ii_tPRP).description=[group_label{groupNo} ' Post  ' week_group_label{week_group}]; 
    end
    x_pos=x_pos+1;
end
xticks([0 1 2 3 5 6 7 8 10 11 12 13]);
xticklabels({'Pre', 'Post', 'Pre', 'Post', 'Pre', 'Post', 'Pre', 'Post',...
    'Pre', 'Post', 'Pre', 'Post' });
ylim([-6 12])
text(-0.5, 10, 'Exo Weeks 1&2', 'FontSize', 14, 'FontWeight', 'bold')
text(5, 10, 'Exo Week 3', 'FontSize', 14, 'FontWeight', 'bold')
text(9.5, 10, 'Poly dA:dT Week 4', 'FontSize', 14, 'FontWeight', 'bold')
text(11,-2,'Mock', 'FontSize', 14, 'FontWeight', 'bold', 'Color',our_colors(1).color)
text(11,-3,'VZV', 'FontSize', 14, 'FontWeight', 'bold', 'Color',our_colors(2).color)

%Now do the glm
fprintf(1, ['\n\nglm for trough delta tPRP olfactory bulb\n'])
tbl = table(glm_tPRP.data',glm_tPRP.group',glm_tPRP.pre_post',glm_tPRP.week',...
    'VariableNames',{'tPRP','group','pre_post','week'});
mdl = fitglm(tbl,'tPRP~group+pre_post+week+group*pre_post*week'...
    ,'CategoricalVars',[2,3,4])


fprintf(1, ['\n\nRanksum or t-test for trough delta tPRP olfactory bulb \n'])
 
%--------------------------------------------------------------------------
% Graph 1 & 2: Delta Peak/Trough tPRP by Trial (Pre/Post Weeks)
%--------------------------------------------------------------------------

plot_delta_tPRP_by_trial = 1; % Enable/disable this section

if plot_delta_tPRP_by_trial == 1
    for ii_electrode = 1:4 % Loop through electrodes
        figNo = figNo + 1;
        hFig = figure(figNo);
        set(hFig, 'units', 'normalized', 'position', [.1 .1 .8 .4]); % Adjust figure size

        subplot(1, 2, 1); hold on; % Peak
        title(['Delta Peak tPRP by Trial - ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}]);
        ylabel('Delta tPRP (dB)');
        % xticks(1:length(trials));
        % xticklabels(trials);
        % xtickangle(45);

        subplot(1, 2, 2); hold on; % Trough
        title(['Delta Trough tPRP by Trial - ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}]);
        ylabel('Delta tPRP (dB)');
        % xticks(1:length(trials));
        % xticklabels(trials);
        % xtickangle(45);

        for trial_idx = 1:length(trials)
            trial_name = trials{trial_idx};

            % Extract week number from trial name (e.g., 'pre-week1' -> 1)
            week_number = str2double(regexp(trial_name, '\d+', 'match'));

            for groupNo = 1:length(groups)
                % Collect data for this group and trial
                these_delta_peak_tPRP = [];
                these_delta_trough_tPRP = [];

                for fileNo = 1:handles.no_files
                    if (handles.group_no(fileNo) == groupNo)&(handles.session_no(fileNo) == trial_idx)

                        %  %Check if this file contains the appropriate week number
                        %  file_week_number = handles.file(fileNo).week_number; %Assumes you have a week_number field in your handles.file struct
                        % 
                        % %Skip files from the wrong week
                        %  if file_week_number ~= week_number
                        %      continue;
                        %  end

                        this_filename = handles.FileName{fileNo};
                        FileName_out = ['OBtoHIPtPRP_' this_filename(11:end)];

                        if exist([PathName_out FileName_out], 'file') %Take out
                            load([PathName_out FileName_out]);

                            add_electrode = isempty(handles.file(fileNo).exclude_electrodes) || ...
                                            ~ismember(actual_electrode(ii_electrode), handles.file(fileNo).exclude_electrodes);

                            if add_electrode
                                log_P_timecourse = handles_out.electrode(ii_electrode).handles_out.mean_peakPower;
                                time = handles_out.electrode(ii_electrode).handles_out.time;

                                % Calculate delta tPRP (Laser - Pre)
                                delta_peak_tPRP = mean(log_P_timecourse((time >= t_ref_end) & (time <= t_laser_end))) - ...
                                                  mean(log_P_timecourse(time <= t_ref_end));

                                delta_trough_tPRP =  mean(handles_out.electrode(ii_electrode).handles_out.mean_troughPower((time >= t_ref_end) & (time <= t_laser_end))) - ...
                                                  mean(handles_out.electrode(ii_electrode).handles_out.mean_troughPower(time <= t_ref_end));


                                these_delta_peak_tPRP = [these_delta_peak_tPRP, delta_peak_tPRP];
                                these_delta_trough_tPRP = [these_delta_trough_tPRP, delta_trough_tPRP];
                            end
                        end
                    end
                end

                % Plotting
               
                x_pos = 7*rem(trial_idx-1,4) + (groupNo - 1)*3 + int8((trial_idx-1)/4);

                if ~isempty(these_delta_peak_tPRP)
                    subplot(1, 2, 1);
                    bar(x_pos, mean(these_delta_peak_tPRP), 0.7, 'FaceColor', our_colors(groupNo).color, 'EdgeColor', 'none'); % Use group color
                    errorbar(x_pos, mean(these_delta_peak_tPRP), std(these_delta_peak_tPRP)/sqrt(length(these_delta_peak_tPRP)), 'k', 'LineStyle', 'none'); % Error bars
                end

                 if ~isempty(these_delta_trough_tPRP)
                    subplot(1, 2, 2);
                    bar(x_pos, mean(these_delta_trough_tPRP), 0.7, 'FaceColor', our_colors(groupNo).color, 'EdgeColor', 'none'); % Use group color
                    errorbar(x_pos, mean(these_delta_trough_tPRP), std(these_delta_trough_tPRP)/sqrt(length(these_delta_trough_tPRP)), 'k', 'LineStyle', 'none'); % Error bars
                 end

                 pfft=1;
            end
        end

        subplot(1, 2, 1);
        xticks([0 1 3 4 7 8 10 11 15 16 18 19 22 23 25 26]);
        xticklabels({'pre-week 1', 'pos-week 1', 'pre-week 1', 'pos-week 1', 'pre-week 2', 'pos-week 2', 'pre-week 2', 'pos-week 2',...
           'pre-week 3', 'pos-week 3', 'pre-week 3', 'pos-week 3', 'pre-week 4', 'pos-week 4', 'pre-week 4', 'pos-week 4' });
        xtickangle(45);
        text(22,3.2, 'mock','Color', '#FFA500');
        text(22,3.8,'vzv', 'Color', '#87CEEB');

        subplot(1, 2, 2);
        xticks([0 1 3 4 7 8 10 11 15 16 18 19 22 23 25 26]);
        xticklabels({'pre-week 1', 'pos-week 1', 'pre-week 1', 'pos-week 1', 'pre-week 2', 'pos-week 2', 'pre-week 2', 'pos-week 2',...
           'pre-week 3', 'pos-week 3', 'pre-week 3', 'pos-week 3', 'pre-week 4', 'pos-week 4', 'pre-week 4', 'pos-week 4' });
        xtickangle(45);
        text(22,3.2, 'mock','Color', '#FFA500');
        text(22,3.8,'vzv', 'Color', '#87CEEB');
    end
end


%--------------------------------------------------------------------------
% Graph 3 & 4: Delta Peak/Trough tPRP Pre-Week 1 vs. Post-Week 4
%--------------------------------------------------------------------------
% 
% plot_delta_tPRP_pre_post = 1; % Enable/disable this section
% 
% if plot_delta_tPRP_pre_post == 1
%     for ii_electrode = 1:4
%         figNo = figNo + 1;
%         hFig = figure(figNo);
%         set(hFig, 'units', 'normalized', 'position', [.1 .1 .4 .4]);
% 
%         subplot(1, 2, 1); hold on; % Peak
%         title(['Delta Peak tPRP Pre-Week 1 vs. Post-Week 4 - ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}]);
%         ylabel('Delta tPRP (dB)');
%         xticks([1 2]);
%         xticklabels({'Pre-Week 1', 'Post-Week 4'});
% 
%         subplot(1, 2, 2); hold on; % Trough
%         title(['Delta Trough tPRP Pre-Week 1 vs. Post-Week 4 - ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}]);
%         ylabel('Delta tPRP (dB)');
%         xticks([1 2]);
%         xticklabels({'Pre-Week 1', 'Post-Week 4'});
% 
% 
%         for groupNo = 1:length(groups)
%             these_delta_peak_tPRP_pre = [];
%             these_delta_peak_tPRP_post = [];
%             these_delta_trough_tPRP_pre = [];
%             these_delta_trough_tPRP_post = [];
% 
% 
%             for fileNo = 1:handles.no_files
%                 if handles.group_no(fileNo) == groupNo
% 
%                     this_filename = handles.FileName{fileNo};
%                     FileName_out = ['OBtoHIPtPRP_' this_filename(11:end)];
%                      file_week_number = handles.file(fileNo).week_number; %Assumes you have a week_number field in your handles.file struct
% 
%                     if exist([PathName_out FileName_out], 'file')
%                         load([PathName_out FileName_out]);
% 
%                         add_electrode = isempty(handles.file(fileNo).exclude_electrodes) || ...
%                                         ~ismember(actual_electrode(ii_electrode), handles.file(fileNo).exclude_electrodes);
% 
%                         if add_electrode
%                             log_P_timecourse = handles_out.electrode(ii_electrode).handles_out.mean_peakPower;
%                             trough_P_timecourse = handles_out.electrode(ii_electrode).handles_out.mean_troughPower;
%                             time = handles_out.electrode(ii_electrode).handles_out.time;
% 
%                             %Check if this file is Pre-Week 1
%                             if file_week_number == 1
% 
%                                 delta_peak_tPRP = mean(log_P_timecourse((time >= t_ref_end) & (time <= t_laser_end))) - ...
%                                                   mean(log_P_timecourse(time <= t_ref_end));
%                                 delta_trough_tPRP = mean(trough_P_timecourse((time >= t_ref_end) & (time <= t_laser_end))) - ...
%                                                   mean(trough_P_timecourse(time <= t_ref_end));
% 
% 
%                                 these_delta_peak_tPRP_pre = [these_delta_peak_tPRP_pre, delta_peak_tPRP];
%                                 these_delta_trough_tPRP_pre = [these_delta_trough_tPRP_pre, delta_trough_tPRP];
%                             end
%                              %Check if this file is Post-Week 4
%                             if file_week_number == 8
% 
%                                 delta_peak_tPRP = mean(log_P_timecourse((time >= t_ref_end) & (time <= t_laser_end))) - ...
%                                                   mean(log_P_timecourse(time <= t_ref_end));
%                                 delta_trough_tPRP = mean(trough_P_timecourse((time >= t_ref_end) & (time <= t_laser_end))) - ...
%                                                   mean(trough_P_timecourse(time <= t_ref_end));
% 
% 
%                                 these_delta_peak_tPRP_post = [these_delta_peak_tPRP_post, delta_peak_tPRP];
%                                 these_delta_trough_tPRP_post = [these_delta_trough_tPRP_post, delta_trough_tPRP];
%                             end
%                         end
%                     end
%                 end
%             end
% 
%             %Plotting
%             % Offset for groups
%             if groupNo == 1
%                 group_offset = -0.2;
%             else
%                 group_offset = 0.2;
%             end
% 
% 
%             %Pre-Week 1
%             x_pos_pre = 1 + group_offset;
%             if ~isempty(these_delta_peak_tPRP_pre)
%                 subplot(1, 2, 1);
%                 bar(x_pos_pre, mean(these_delta_peak_tPRP_pre), 0.3, 'FaceColor', our_colors(groupNo).color, 'EdgeColor', 'none');
%                 errorbar(x_pos_pre, mean(these_delta_peak_tPRP_pre), std(these_delta_peak_tPRP_pre)/sqrt(length(these_delta_peak_tPRP_pre)), 'k', 'LineStyle', 'none');
%             end
%             if ~isempty(these_delta_trough_tPRP_pre)
%                 subplot(1, 2, 2);
%                 bar(x_pos_pre, mean(these_delta_trough_tPRP_pre), 0.3, 'FaceColor', our_colors(groupNo).color, 'EdgeColor', 'none');
%                 errorbar(x_pos_pre, mean(these_delta_trough_tPRP_pre), std(these_delta_trough_tPRP_pre)/sqrt(length(these_delta_trough_tPRP_pre)), 'k', 'LineStyle', 'none');
%             end
%             %Post-Week 4
%             x_pos_post = 2 + group_offset;
%              if ~isempty(these_delta_peak_tPRP_post)
%                 subplot(1, 2, 1);
%                 bar(x_pos_post, mean(these_delta_peak_tPRP_post), 0.3, 'FaceColor', our_colors(groupNo).color, 'EdgeColor', 'none');
%                 errorbar(x_pos_post, mean(these_delta_peak_tPRP_post), std(these_delta_peak_tPRP_post)/sqrt(length(these_delta_peak_tPRP_post)), 'k', 'LineStyle', 'none');
%             end
%             if ~isempty(these_delta_trough_tPRP_post)
%                 subplot(1, 2, 2);
%                 bar(x_pos_post, mean(these_delta_trough_tPRP_post), 0.3, 'FaceColor', our_colors(groupNo).color, 'EdgeColor', 'none');
%                 errorbar(x_pos_post, mean(these_delta_trough_tPRP_post), std(these_delta_trough_tPRP_post)/sqrt(length(these_delta_trough_tPRP_post)), 'k', 'LineStyle', 'none');
%             end
%         end
% 
%         subplot(1, 2, 1);
%         legend(group_label{1:length(groups)});
% 
%         subplot(1, 2, 2);
%         legend(group_label{1:length(groups)});
% 
%     end
% end

%--------------------------------------------------------------------------
% Graph 5 & 6: tPRP Peak/Trough with Pre/Laser Windows by Trial
%--------------------------------------------------------------------------

if these_plots_shown(5) == 1 
    for ii_electrode = 1:4
            figNo = figNo + 1;
            hFigPeak = figure(figNo);
            set(hFigPeak, 'units', 'normalized', 'position', [.1 .1 .8 .4]);
            subplot(1, 2, 1);
            title(['tPRP Peak by Trial pre Laser- ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}]);
            hold on
            subplot(1, 2,2 );
            title(['tPRP Trough by Trial pre Laser- ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}]);
            hold on
       
        for trial_idx = 1:length(trials)
            trial_name = trials{trial_idx};

           
            for groupNo = 1:length(groups)
                % Collect data for this group and trial
                these_peak_tPRP_pre = [];
                these_peak_tPRP_laser = [];
                these_trough_tPRP_pre = [];
                these_trough_tPRP_laser = [];


                for fileNo = 1:handles.no_files
                      if (handles.group_no(fileNo) == groupNo)&(handles.session_no(fileNo) == trial_idx)
                        this_filename = handles.FileName{fileNo};
                         file_week_number = handles.file(fileNo).week_number; %Assumes you have a week_number field in your handles.file struct

                        FileName_out = ['OBtoHIPtPRP_' this_filename(11:end)];

                        if exist([PathName_out FileName_out], 'file')

                     
                            load([PathName_out FileName_out]);

                            add_electrode = isempty(handles.file(fileNo).exclude_electrodes) || ...
                                            ~ismember(actual_electrode(ii_electrode), handles.file(fileNo).exclude_electrodes);

                            if add_electrode
                                log_P_timecourse = handles_out.electrode(ii_electrode).handles_out.mean_peakPower;
                                trough_P_timecourse = handles_out.electrode(ii_electrode).handles_out.mean_troughPower;
                                time = handles_out.electrode(ii_electrode).handles_out.time;

                                % Extract tPRP for pre-laser and laser windows
                                peak_tPRP_pre = mean(log_P_timecourse(time <= t_ref_end));
                                peak_tPRP_laser = mean(log_P_timecourse((time >= t_ref_end) & (time <= t_laser_end)));
                                trough_tPRP_pre = mean(trough_P_timecourse(time <= t_ref_end));
                                trough_tPRP_laser = mean(trough_P_timecourse((time >= t_ref_end) & (time <= t_laser_end)));


                                these_peak_tPRP_pre = [these_peak_tPRP_pre, peak_tPRP_pre];
                                these_peak_tPRP_laser = [these_peak_tPRP_laser, peak_tPRP_laser];
                                these_trough_tPRP_pre = [these_trough_tPRP_pre, trough_tPRP_pre];
                                these_trough_tPRP_laser = [these_trough_tPRP_laser, trough_tPRP_laser];
                            end
                        end
                    end
                end

                % Plotting
                % Offset for groups
                if groupNo == 1
                    group_offset = -0.15;
                else
                    group_offset = 0.15;
                end
                  x_pos = 7*rem(trial_idx-1,4) + (groupNo - 1)*3 + int8((trial_idx-1)/4);
                  if ~isempty(these_delta_peak_tPRP)
                    subplot(1, 2, 1);
                    bar(x_pos, mean(these_peak_tPRP_pre), 0.7, 'FaceColor', our_colors(groupNo).color, 'EdgeColor', 'none'); % Use group color
                    errorbar(x_pos, mean(these_peak_tPRP_pre), std(these_peak_tPRP_pre)/sqrt(length(these_peak_tPRP_pre)), 'k', 'LineStyle', 'none'); % Error bars
                end

                 if ~isempty(these_delta_trough_tPRP)
                    subplot(1, 2, 2);
                    bar(x_pos, mean(these_trough_tPRP_pre), 0.7, 'FaceColor', our_colors(groupNo).color, 'EdgeColor', 'none'); % Use group color
                    errorbar(x_pos, mean(these_trough_tPRP_pre), std(these_trough_tPRP_pre)/sqrt(length(these_trough_tPRP_pre)), 'k', 'LineStyle', 'none'); % Error bars
                 end
                 % if these_plots_shown(5) == 1
                 %    % Peak tPRP
                 %    subplot(2, length(trials)/2, trial_idx, 'Parent', hFigPeak); hold on; % Adjust subplot layout as needed
                 %    bar(x_pos, mean(these_peak_tPRP_pre), 0.15, 'FaceColor', our_colors(groupNo).color * 0.7, 'EdgeColor', 'none'); % Pre - lighter color
                 %    % bar(x_pos + 0.1, mean(these_peak_tPRP_laser), 0.15, 'FaceColor', our_colors(groupNo).color, 'EdgeColor', 'none'); % Laser - original color
                 %    errorbar(x_pos, mean(these_peak_tPRP_pre), std(these_peak_tPRP_pre)/sqrt(length(these_peak_tPRP_pre)), 'k', 'LineStyle', 'none');
                 %    % errorbar(x_pos + 0.1, mean(these_peak_tPRP_laser), std(these_peak_tPRP_laser)/sqrt(length(these_peak_tPRP_laser)), 'k', 'LineStyle', 'none');
                 %    % xticks(trial_idx);
                 %    % xticklabels(trial_name);
                 %    % xtickangle(45);
                 %    ylabel('tPRP (dB)');
                 % 
                 % end
                 % 
                 % if these_plots_shown(7) == 1
                 %    % Trough tPRP
                 %    subplot(2, length(trials)/2, trial_idx, 'Parent', hFigTrough); hold on; % Adjust subplot layout as needed
                 %    bar(x_pos, mean(these_trough_tPRP_pre), 0.15, 'FaceColor', our_colors(groupNo).color * 0.7, 'EdgeColor', 'none'); % Pre - lighter color
                 %    %bar(x_pos + 0.1, mean(these_trough_tPRP_laser), 0.15, 'FaceColor', our_colors(groupNo).color, 'EdgeColor', 'none'); % Laser - original color
                 %    errorbar(x_pos, mean(these_trough_tPRP_pre), std(these_trough_tPRP_pre)/sqrt(length(these_trough_tPRP_pre)), 'k', 'LineStyle', 'none');
                 %    %errorbar(x_pos + 0.1, mean(these_trough_tPRP_laser), std(these_trough_tPRP_laser)/sqrt(length(these_trough_tPRP_laser)), 'k', 'LineStyle', 'none');
                 %    % xticks(trial_idx);
                 %    % xticklabels(trial_name);
                 %    % xtickangle(45);
                 %    ylabel('tPRP (dB)');
                 %    % title(trial_name);
                 % end
                 pff=1;
            end

             if these_plots_shown(5) == 1
                % subplot(2, length(trials)/2, 1, 'Parent', hFigPeak); % Put legend on the first subplot
                subplot(1, 2, 1);
                xticks([0 1 3 4 7 8 10 11 15 16 18 19 22 23 25 26]);
                xticklabels({'pre-week 1', 'pos-week 1', 'pre-week 1', 'pos-week 1', 'pre-week 2', 'pos-week 2', 'pre-week 2', 'pos-week 2',...
                    'pre-week 3', 'pos-week 3', 'pre-week 3', 'pos-week 3', 'pre-week 4', 'pos-week 4', 'pre-week 4', 'pos-week 4' });
                xtickangle(45);
                text(22,56, 'mock','Color', '#FFA500');
                text(22,58,'vzv', 'Color', '#87CEEB');

                subplot(1, 2, 2);
                xticks([0 1 3 4 7 8 10 11 15 16 18 19 22 23 25 26]);
                xticklabels({'pre-week 1', 'pos-week 1', 'pre-week 1', 'pos-week 1', 'pre-week 2', 'pos-week 2', 'pre-week 2', 'pos-week 2',...
                    'pre-week 3', 'pos-week 3', 'pre-week 3', 'pos-week 3', 'pre-week 4', 'pos-week 4', 'pre-week 4', 'pos-week 4' });
                xtickangle(45);
                text(22,56, 'mock','Color', '#FFA500');
                text(22,58,'vzv', 'Color', '#87CEEB');

             end

             % if these_plots_shown(7) == 1
             %    subplot(2, length(trials)/2, 1, 'Parent', hFigTrough); % Put legend on the first subplot
             %    legend(group_label{1:length(groups)});
             % end

        end
    end
end


if these_plots_shown(5) == 1 
    for ii_electrode = 1:4
            figNo = figNo + 1;
            hFigPeak = figure(figNo);
            set(hFigPeak, 'units', 'normalized', 'position', [.1 .1 .8 .4]);
            subplot(1, 2, 1);
            title(['tPRP Peak by Trial Laser- ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}]);
            hold on
            subplot(1, 2,2 );
            title(['tPRP Trough by Trial Laser- ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}]);
            hold on
       
        for trial_idx = 1:length(trials)
            trial_name = trials{trial_idx};

           
            for groupNo = 1:length(groups)
                % Collect data for this group and trial
                these_peak_tPRP_pre = [];
                these_peak_tPRP_laser = [];
                these_trough_tPRP_pre = [];
                these_trough_tPRP_laser = [];


                for fileNo = 1:handles.no_files
                      if (handles.group_no(fileNo) == groupNo)&(handles.session_no(fileNo) == trial_idx)
                        this_filename = handles.FileName{fileNo};
                         file_week_number = handles.file(fileNo).week_number; %Assumes you have a week_number field in your handles.file struct

                        FileName_out = ['OBtoHIPtPRP_' this_filename(11:end)];

                        if exist([PathName_out FileName_out], 'file')

                     
                            load([PathName_out FileName_out]);

                            add_electrode = isempty(handles.file(fileNo).exclude_electrodes) || ...
                                            ~ismember(actual_electrode(ii_electrode), handles.file(fileNo).exclude_electrodes);

                            if add_electrode
                                log_P_timecourse = handles_out.electrode(ii_electrode).handles_out.mean_peakPower;
                                trough_P_timecourse = handles_out.electrode(ii_electrode).handles_out.mean_troughPower;
                                time = handles_out.electrode(ii_electrode).handles_out.time;

                                % Extract tPRP for pre-laser and laser windows
                                peak_tPRP_pre = mean(log_P_timecourse(time <= t_ref_end));
                                peak_tPRP_laser = mean(log_P_timecourse((time >= t_ref_end) & (time <= t_laser_end)));
                                trough_tPRP_pre = mean(trough_P_timecourse(time <= t_ref_end));
                                trough_tPRP_laser = mean(trough_P_timecourse((time >= t_ref_end) & (time <= t_laser_end)));


                                these_peak_tPRP_pre = [these_peak_tPRP_pre, peak_tPRP_pre];
                                these_peak_tPRP_laser = [these_peak_tPRP_laser, peak_tPRP_laser];
                                these_trough_tPRP_pre = [these_trough_tPRP_pre, trough_tPRP_pre];
                                these_trough_tPRP_laser = [these_trough_tPRP_laser, trough_tPRP_laser];
                            end
                        end
                    end
                end

                % Plotting
                % Offset for groups
                if groupNo == 1
                    group_offset = -0.15;
                else
                    group_offset = 0.15;
                end
                  x_pos = 7*rem(trial_idx-1,4) + (groupNo - 1)*3 + int8((trial_idx-1)/4);
                  if ~isempty(these_delta_peak_tPRP)
                    subplot(1, 2, 1);
                    bar(x_pos, mean(these_peak_tPRP_laser), 0.7, 'FaceColor', our_colors(groupNo).color, 'EdgeColor', 'none'); % Use group color
                    errorbar(x_pos, mean(these_peak_tPRP_laser), std(these_peak_tPRP_laser)/sqrt(length(these_peak_tPRP_laser)), 'k', 'LineStyle', 'none'); % Error bars
                end

                 if ~isempty(these_delta_trough_tPRP)
                    subplot(1, 2, 2);
                    bar(x_pos, mean(these_trough_tPRP_laser), 0.7, 'FaceColor', our_colors(groupNo).color, 'EdgeColor', 'none'); % Use group color
                    errorbar(x_pos, mean(these_trough_tPRP_laser), std(these_trough_tPRP_laser)/sqrt(length(these_trough_tPRP_laser)), 'k', 'LineStyle', 'none'); % Error bars
                 end
                 % if these_plots_shown(5) == 1
                 %    % Peak tPRP
                 %    subplot(2, length(trials)/2, trial_idx, 'Parent', hFigPeak); hold on; % Adjust subplot layout as needed
                 %    bar(x_pos, mean(these_peak_tPRP_pre), 0.15, 'FaceColor', our_colors(groupNo).color * 0.7, 'EdgeColor', 'none'); % Pre - lighter color
                 %    % bar(x_pos + 0.1, mean(these_peak_tPRP_laser), 0.15, 'FaceColor', our_colors(groupNo).color, 'EdgeColor', 'none'); % Laser - original color
                 %    errorbar(x_pos, mean(these_peak_tPRP_pre), std(these_peak_tPRP_pre)/sqrt(length(these_peak_tPRP_pre)), 'k', 'LineStyle', 'none');
                 %    % errorbar(x_pos + 0.1, mean(these_peak_tPRP_laser), std(these_peak_tPRP_laser)/sqrt(length(these_peak_tPRP_laser)), 'k', 'LineStyle', 'none');
                 %    % xticks(trial_idx);
                 %    % xticklabels(trial_name);
                 %    % xtickangle(45);
                 %    ylabel('tPRP (dB)');
                 % 
                 % end
                 % 
                 % if these_plots_shown(7) == 1
                 %    % Trough tPRP
                 %    subplot(2, length(trials)/2, trial_idx, 'Parent', hFigTrough); hold on; % Adjust subplot layout as needed
                 %    bar(x_pos, mean(these_trough_tPRP_pre), 0.15, 'FaceColor', our_colors(groupNo).color * 0.7, 'EdgeColor', 'none'); % Pre - lighter color
                 %    %bar(x_pos + 0.1, mean(these_trough_tPRP_laser), 0.15, 'FaceColor', our_colors(groupNo).color, 'EdgeColor', 'none'); % Laser - original color
                 %    errorbar(x_pos, mean(these_trough_tPRP_pre), std(these_trough_tPRP_pre)/sqrt(length(these_trough_tPRP_pre)), 'k', 'LineStyle', 'none');
                 %    %errorbar(x_pos + 0.1, mean(these_trough_tPRP_laser), std(these_trough_tPRP_laser)/sqrt(length(these_trough_tPRP_laser)), 'k', 'LineStyle', 'none');
                 %    % xticks(trial_idx);
                 %    % xticklabels(trial_name);
                 %    % xtickangle(45);
                 %    ylabel('tPRP (dB)');
                 %    % title(trial_name);
                 % end
                 pff=1;
            end

             if these_plots_shown(5) == 1
                % subplot(2, length(trials)/2, 1, 'Parent', hFigPeak); % Put legend on the first subplot
                subplot(1, 2, 1);
                xticks([0 1 3 4 7 8 10 11 15 16 18 19 22 23 25 26]);
                xticklabels({'pre-week 1', 'pos-week 1', 'pre-week 1', 'pos-week 1', 'pre-week 2', 'pos-week 2', 'pre-week 2', 'pos-week 2',...
                    'pre-week 3', 'pos-week 3', 'pre-week 3', 'pos-week 3', 'pre-week 4', 'pos-week 4', 'pre-week 4', 'pos-week 4' });
                xtickangle(45);
                text(22,66, 'mock','Color', '#FFA500');
                text(22,68,'vzv', 'Color', '#87CEEB');

                subplot(1, 2, 2);
                xticks([0 1 3 4 7 8 10 11 15 16 18 19 22 23 25 26]);
                xticklabels({'pre-week 1', 'pos-week 1', 'pre-week 1', 'pos-week 1', 'pre-week 2', 'pos-week 2', 'pre-week 2', 'pos-week 2',...
                    'pre-week 3', 'pos-week 3', 'pre-week 3', 'pos-week 3', 'pre-week 4', 'pos-week 4', 'pre-week 4', 'pos-week 4' });
                xtickangle(45);
               text(22,66, 'mock','Color', '#FFA500');
               text(22,68,'vzv', 'Color', '#87CEEB');
               

             end
        end
    end
end

%--------------------------------------------------------------------------
% Graph 7 & 8: tPRP Peak/Trough with Pre/Laser Windows Pre-Week 1 vs. Post-Week 4
%--------------------------------------------------------------------------
% 
% plot_tPRP_pre_post_windows = 1; % Enable/disable this section
% 
% if plot_tPRP_pre_post_windows == 1
%     for ii_electrode = 1:4
%         figNo = figNo + 1;
%         hFig = figure(figNo);
%         set(hFig, 'units', 'normalized', 'position', [.1 .1 .6 .4]); % Adjust figure size for two subplots
% 
%         subplot(1, 2, 1); hold on; % Peak
%         title(['tPRP Peak Pre/Laser - Pre-Week 1 vs. Post-Week 4 - ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}]);
%         ylabel('tPRP (dB)');
%         xticks([1 2]);
%         xticklabels({'Pre-Week 1', 'Post-Week 4'});
% 
%         subplot(1, 2, 2); hold on; % Trough
%         title(['tPRP Trough Pre/Laser - Pre-Week 1 vs. Post-Week 4 - ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}]);
%         ylabel('tPRP (dB)');
%         xticks([1 2]);
%         xticklabels({'Pre-Week 1', 'Post-Week 4'});
% 
%         for groupNo = 1:length(groups)
%             these_peak_tPRP_pre_pre = [];
%             these_peak_tPRP_laser_pre = [];
%             these_peak_tPRP_pre_post = [];
%             these_peak_tPRP_laser_post = [];
%             these_trough_tPRP_pre_pre = [];
%             these_trough_tPRP_laser_pre = [];
%             these_trough_tPRP_pre_post = [];
%             these_trough_tPRP_laser_post = [];
% 
% 
%             for fileNo = 1:handles.no_files
%                 if handles.group_no(fileNo) == groupNo
%                     this_filename = handles.FileName{fileNo};
%                     FileName_out = ['OBtoHIPtPRP_' this_filename(11:end)];
%                     file_week_number = handles.file(fileNo).week_number; %Assumes you have a week_number field in your handles.file struct
% 
%                     if exist([PathName_out FileName_out], 'file')
%                         load([PathName_out FileName_out]);
% 
%                         add_electrode = isempty(handles.file(fileNo).exclude_electrodes) || ...
%                                         ~ismember(actual_electrode(ii_electrode), handles.file(fileNo).exclude_electrodes);
% 
%                         if add_electrode
%                             log_P_timecourse = handles_out.electrode(ii_electrode).handles_out.mean_peakPower;
%                             trough_P_timecourse = handles_out.electrode(ii_electrode).handles_out.mean_troughPower;
%                             time = handles_out.electrode(ii_electrode).handles_out.time;
% 
%                             %Check if this file is Pre-Week 1
%                             if file_week_number == 1
%                                 peak_tPRP_pre = mean(log_P_timecourse(time <= t_ref_end));
%                                 peak_tPRP_laser = mean(log_P_timecourse((time >= t_ref_end) & (time <= t_laser_end)));
%                                 trough_tPRP_pre = mean(trough_P_timecourse(time <= t_ref_end));
%                                 trough_tPRP_laser = mean(trough_P_timecourse((time >= t_ref_end) & (time <= t_laser_end)));
% 
%                                 these_peak_tPRP_pre_pre = [these_peak_tPRP_pre_pre, peak_tPRP_pre];
%                                 these_peak_tPRP_laser_pre = [these_peak_tPRP_laser_pre, peak_tPRP_laser];                               
%                                 these_trough_tPRP_pre_pre = [these_trough_tPRP_pre_pre, trough_tPRP_pre];
%                                 these_trough_tPRP_laser_pre = [these_trough_tPRP_laser_pre, trough_tPRP_laser];
%                             end
%                              %Check if this file is Post-Week 4
%                             if file_week_number == 8
%                                 peak_tPRP_pre = mean(log_P_timecourse(time <= t_ref_end));
%                                 peak_tPRP_laser = mean(log_P_timecourse((time >= t_ref_end) & (time <= t_laser_end)));
%                                 trough_tPRP_pre = mean(trough_P_timecourse(time <= t_ref_end));
%                                 trough_tPRP_laser = mean(trough_P_timecourse((time >= t_ref_end) & (time <= t_laser_end)));
% 
% 
%                                 these_peak_tPRP_pre_post = [these_peak_tPRP_pre_post, peak_tPRP_pre];
%                                 these_peak_tPRP_laser_post = [these_peak_tPRP_laser_post, peak_tPRP_laser];
%                                 these_trough_tPRP_pre_post = [these_trough_tPRP_pre_post, trough_tPRP_pre];
%                                 these_trough_tPRP_laser_post = [these_trough_tPRP_laser_post, trough_tPRP_laser];
%                             end
%                         end
%                     end
%                 end
%             end
% 
%             %Plotting
%             % Offset for groups
%             if groupNo == 1
%                 group_offset = -0.2;
%             else
%                 group_offset = 0.2;
%             end
% 
%             %Pre-Week 1
%             x_pos_pre = 1 + group_offset;
%             subplot(1, 2, 1);
%             hold on;
%             if ~isempty(these_peak_tPRP_pre_pre)
%                 bar(x_pos_pre - 0.1, mean(these_peak_tPRP_pre_pre), 0.15, 'FaceColor', our_colors(groupNo).color * 0.7, 'EdgeColor', 'none'); % Pre - lighter color
%                 errorbar(x_pos_pre - 0.1, mean(these_peak_tPRP_pre_pre), std(these_peak_tPRP_pre_pre)/sqrt(length(these_peak_tPRP_pre_pre)), 'k', 'LineStyle', 'none');
%             end
%             if ~isempty(these_peak_tPRP_laser_pre)
%                 bar(x_pos_pre + 0.1, mean(these_peak_tPRP_laser_pre), 0.15, 'FaceColor', our_colors(groupNo).color, 'EdgeColor', 'none'); % Laser - original color
%                 errorbar(x_pos_pre + 0.1, mean(these_peak_tPRP_laser_pre), std(these_peak_tPRP_laser_pre)/sqrt(length(these_peak_tPRP_laser_pre)), 'k', 'LineStyle', 'none');
%             end
% 
%             subplot(1, 2, 2);
%             hold on;
%             if ~isempty(these_trough_tPRP_pre_pre)
%                 bar(x_pos_pre - 0.1, mean(these_trough_tPRP_pre_pre), 0.15, 'FaceColor', mean(these_trough_tPRP_pre_pre), 0.15, 'FaceColor', our_colors(groupNo).color * 0.7, 'EdgeColor', 'none'); % Pre - lighter color
%                 errorbar(x_pos_pre - 0.1, mean(these_trough_tPRP_pre_pre), std(these_trough_tPRP_pre_pre)/sqrt(length(these_trough_tPRP_pre_pre)), 'k', 'LineStyle', 'none');
%             end
%             if ~isempty(these_trough_tPRP_laser_pre)
%                 bar(x_pos_pre + 0.1, mean(these_trough_tPRP_laser_pre), 0.15, 'FaceColor', our_colors(groupNo).color, 'EdgeColor', 'none'); % Laser - original color
%                 errorbar(x_pos_pre + 0.1, mean(these_trough_tPRP_laser_pre), std(these_trough_tPRP_laser_pre)/sqrt(length(these_trough_tPRP_laser_pre)), 'k', 'LineStyle', 'none');
%             end
% 
%             %Post-Week 4
%             x_pos_post = 2 + group_offset;
%             subplot(1, 2, 1);
%             if ~isempty(these_peak_tPRP_pre_post)
%                 bar(x_pos_post - 0.1, mean(these_peak_tPRP_pre_post), 0.15, 'FaceColor', our_colors(groupNo).color * 0.7, 'EdgeColor', 'none'); % Pre - lighter color
%                 errorbar(x_pos_post - 0.1, mean(these_peak_tPRP_pre_post), std(these_peak_tPRP_pre_post)/sqrt(length(these_peak_tPRP_pre_post)), 'k', 'LineStyle', 'none');
%             end
%             if ~isempty(these_peak_tPRP_laser_post)
%                 bar(x_pos_post + 0.1, mean(these_peak_tPRP_laser_post), 0.15, 'FaceColor', our_colors(groupNo).color, 'EdgeColor', 'none'); % Laser - original color
%                 errorbar(x_pos_post + 0.1, mean(these_peak_tPRP_laser_post), std(these_peak_tPRP_laser_post)/sqrt(length(these_peak_tPRP_laser_post)), 'k', 'LineStyle', 'none');
%             end
% 
%             subplot(1, 2, 2);
%             if ~isempty(these_trough_tPRP_pre_post)
%                 bar(x_pos_post - 0.1, mean(these_trough_tPRP_pre_post), 0.15, 'FaceColor', our_colors(groupNo).color * 0.7, 'EdgeColor', 'none'); % Pre - lighter color
%                 errorbar(x_pos_post - 0.1, mean(these_trough_tPRP_pre_post), std(these_trough_tPRP_pre_post)/sqrt(length(these_trough_tPRP_pre_post)), 'k', 'LineStyle', 'none');
%             end
%             if ~isempty(these_trough_tPRP_laser_post)
%                 bar(x_pos_post + 0.1, mean(these_trough_tPRP_laser_post), 0.15, 'FaceColor', our_colors(groupNo).color, 'EdgeColor', 'none'); % Laser - original color
%                 errorbar(x_pos_post + 0.1, mean(these_trough_tPRP_laser_post), std(these_trough_tPRP_laser_post)/sqrt(length(these_trough_tPRP_laser_post)), 'k', 'LineStyle', 'none');
%             end
%         end
% 
%         subplot(1, 2, 1);
%         legend(group_label{1:length(groups)});
% 
%         subplot(1, 2, 2);
%         legend(group_label{1:length(groups)});
% 
%     end
% end
% 
% %--------------------------------------------------------------------------
% % Remaining Plots (Modulation Index, Angle Variance, etc.)
% %--------------------------------------------------------------------------
% 
% % Add code for plots 9-13 here, modifying as needed
% 
% %Show modulation index bar graphs for each electrode for each time window
% 
% if these_plots_shown(9)==1
% 
% edges=[-30:5:15];
% 
% rand_offset=0.8;
% 
% for ii_electrode=1:4
% 
% figNo=figNo+1;
% 
% try
% 
% close(figNo)
% 
% catch
% 
% end
% 
% hFig2 = figure(figNo);
% 
% set(hFig2, 'units','normalized','position',[.1 .1 .3 .3])
% 
% hold on
% 
% bar_offset=0;
% 
% glm_lfp_power=[];
% 
% glm_ii_lfp=0;
% 
% id_ii=0;
% 
% input_data=[];
% 
% for groupNo=1:length(groups)
% 
% %Find the files for this mouse
% 
% these_files=[];
% 
% for fileNo=1:handles.no_files
% 
% if handles.group_no(fileNo)==groupNo
% 
% %Does this file exist?
% 
% this_filename=handles.FileName{fileNo};
% 
% FileName_out=['OBtoHIPtPRP_' this_filename(11:end)];
% 
% if exist([PathName_out FileName_out])~=0
% 
% these_files=[these_files fileNo];
% 
% end
% 
% end
% 
% end
% 
% %Get MI
% 
% these_delta_log_P_laser=[];
% 
% these_delta_log_P_post=[];
% 
% these_mice=[];
% 
% if ~isempty(these_files)
% 
% for fileNo=these_files
% 
% %Plot the timecourse
% 
% % PathName_out=handles.PathName{fileNo};
% 
% this_filename=handles.FileName{fileNo};
% 
% FileName_out=['OBtoHIPtPRP_' this_filename(11:end)];
% 
% load([PathName_out FileName_out])
% 
% add_electrode=0;
% 
% if isempty(handles.file(fileNo).exclude_electrodes)
% 
% add_electrode=1;
% 
% else
% 
% if sum(handles.file(fileNo).exclude_electrodes==actual_electrode(ii_electrode))==0
% 
% add_electrode=1;
% 
% end
% 
% end
% 
% if add_electrode==1
% 
% modulation_index=handles_out.electrode(ii_electrode).handles_out.modulation_index;
% 
% these_delta_log_P_laser=[these_delta_log_P_laser modulation_index(1)];
% 
% these_delta_log_P_post=[these_delta_log_P_post modulation_index(2)];
% 
% these_mice=[these_mice handles.mouse_no(fileNo)];
% 
% end
% 
% pffft=1;
% 
% end
% 
% these_delta_log_P_laser=these_delta_log_P_laser(~isnan(these_delta_log_P_laser));
% 
% these_mice_laser=these_mice(~isnan(these_delta_log_P_laser));
% 
% these_delta_log_P_post=these_delta_log_P_post(~isnan(these_delta_log_P_post));
% 
% these_mice_post=these_mice(~isnan(these_delta_log_P_post));
% 
% %Laser
% 
% bar(bar_offset,mean(these_delta_log_P_laser),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
% 
% [mean_out, CIout]=drgViolinPoint(these_delta_log_P_laser,edges,bar_offset,rand_offset,'k','k',3);
% 
% glm_lfp_power.data(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=these_delta_log_P_laser;
% 
% glm_lfp_power.time_window(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=0*ones(1,length(these_delta_log_P_laser));
% 
% switch groupNo
% 
% case 1
% 
% glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=1;
% 
% glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=1;
% 
% case 2
% 
% glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=1;
% 
% glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=2;
% 
% case 3
% 
% glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=2;
% 
% glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=1;
% 
% case 4
% 
% glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=2;
% 
% glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=2;
% 
% end
% 
% glm_lfp_power.mouseNo(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=these_mice_laser;
% 
% glm_ii_lfp=glm_ii_lfp+length(these_delta_log_P_laser);
% 
% id_ii=id_ii+1;
% 
% input_data(id_ii).data=these_delta_log_P_laser;
% 
% input_data(id_ii).description=['laser ' group_label{groupNo}];
% 
% %Post
% 
% bar_offset=bar_offset+1;
% 
% bar(bar_offset,mean(these_delta_log_P_post),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
% 
% [mean_out, CIout]=drgViolinPoint(these_delta_log_P_post,edges,bar_offset,rand_offset,'k','k',3);
% 
% glm_lfp_power.data(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=these_delta_log_P_post;
% 
% glm_lfp_power.time_window(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=1*ones(1,length(these_delta_log_P_post));
% 
% switch groupNo
% 
% case 1
% 
% glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=1;
% 
% glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=1;
% 
% case 2
% 
% glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=1;
% 
% glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=2;
% 
% case 3
% 
% glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=2;
% 
% glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=1;
% 
% case 4
% 
% glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=2;
% 
% glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=2;
% 
% end
% 
% glm_lfp_power.mouseNo(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=these_mice_post;
% 
% glm_ii_lfp=glm_ii_lfp+length(these_delta_log_P_post);
% 
% id_ii=id_ii+1;
% 
% input_data(id_ii).data=these_delta_log_P_post;
% 
% input_data(id_ii).description=['post ' group_label{groupNo}];
% 
% bar_offset=bar_offset+2;
% 
% else
% 
% bar_offset=bar_offset+3;
% 
% end
% 
% end
% 
% title(['Modulation Index ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}])
% 
% xticks([0.5 3.5 6.5 9.5])
% 
% xticklabels({group_label{1},group_label{2},group_label{3},group_label{4}})
% 
% % ylim([0 0.03])
% 
% ylabel('Modulation Index')
% 
% %Perform the glm for MI per mouse per odor pair
% 
% fprintf(1, ['glm for Modulation Index for ' handles.electrode_label{handles.peakLFPNo(ii_electrode)} '\n'])
% 
% fprintf(fileID, ['glm for Modulation Index for ' handles.electrode_label{handles.peakLFPNo(ii_electrode)} '\n']);
% 
% % tbl = table(glm_lfp_power.data',glm_lfp_power.time_window',glm_lfp_power.group',...
% 
% % 'VariableNames',{'logP','laser_post','group'});
% 
% % mdl = fitglm(tbl,'logP~laser_post+group+laser_post*group'...
% 
% % ,'CategoricalVars',[2,3])
% 
% % First, create a cell array of character labels corresponding to your group numbers
% 
% % group_labels = {'Control', '5xFAD', 'Group3', 'Group4'}; % Add all your group labels here
% 
% tbl = table(glm_lfp_power.data',categorical(glm_lfp_power.time_window', [0 1], time_window_labels),categorical(glm_lfp_power.genotype',[1], genotype_label),categorical(glm_lfp_power.laser',[1], laser_label),'VariableNames',{'logP','time_window','genotype','laser'});
% 
% mdl = fitglm(tbl,'logP~time_window+genotype+laser+time_window*genotype*laser')
% 
% txt = evalc('mdl');
% 
% txt=regexp(txt,'','split');
% 
% txt=cell2mat(txt);
% 
% txt=regexp(txt,'','split');
% 
% txt=cell2mat(txt);
% 
% fprintf(fileID,'%s\n', txt);
% 
% %Do the ranksum/t-test
% 
% fprintf(1, ['\n\nRanksum or t-test p values for Modulation Index for ' handles.electrode_label{handles.peakLFPNo(ii_electrode)} '\n']);
% 
% fprintf(fileID, ['\n\nRanksum or t-test p values for Modulation Index for ' handles.electrode_label{handles.peakLFPNo(ii_electrode)} '\n']);
% 
% [output_data] = drgMutiRanksumorTtest(input_data,fileID);
% 
% end
% 
% end
% 
% %Show peak angle variance bar graphs for each electrode for each time window
% 
% if these_plots_shown(12)==1
% 
% edges=[-30:5:15];
% 
% rand_offset=0.8;
% 
% for ii_electrode=1:4
% 
% figNo=figNo+1;
% 
% try
% 
% close(figNo)
% 
% catch
% 
% end
% 
% hFig2 = figure(figNo);
% 
% set(hFig2, 'units','normalized','position',[.1 .1 .3 .3])
% 
% hold on
% 
% bar_offset=0;
% 
% glm_lfp_power=[];
% 
% glm_ii_lfp=0;
% 
% id_ii=0;
% 
% input_data=[];
% 
% for groupNo=1:length(groups)
% 
% %Find the files for this mouse
% 
% these_files=[];
% 
% for fileNo=1:handles.no_files
% 
% if handles.group_no(fileNo)==groupNo
% 
% %Does this file exist?
% 
% this_filename=handles.FileName{fileNo};
% 
% FileName_out=['OBtoHIPtPRP_' this_filename(11:end)];
% 
% if exist([PathName_out FileName_out])~=0
% 
% these_files=[these_files fileNo];
% 
% end
% 
% end
% 
% end
% 
% %Get variance
% 
% these_delta_log_P_laser=[];
% 
% these_delta_log_P_post=[];
% 
% these_mice=[];
% 
% if ~isempty(these_files)
% 
% for fileNo=these_files
% 
% %Plot the timecourse
% 
% % PathName_out=handles.PathName{fileNo};
% 
% this_filename=handles.FileName{fileNo};
% 
% FileName_out=['OBtoHIPtPRP_' this_filename(11:end)];
% 
% load([PathName_out FileName_out])
% 
% add_electrode=0;
% 
% if isempty(handles.file(fileNo).exclude_electrodes)
% 
% add_electrode=1;
% 
% else
% 
% if sum(handles.file(fileNo).exclude_electrodes==actual_electrode(ii_electrode))==0
% 
% add_electrode=1;
% 
% end
% 
% end
% 
% if add_electrode==1
% 
% peak_angle_variance=handles_out.electrode(ii_electrode).handles_out.peak_angle_variance;
% 
% these_delta_log_P_laser=[these_delta_log_P_laser peak_angle_variance(1)];
% 
% these_delta_log_P_post=[these_delta_log_P_post peak_angle_variance(2)];
% 
% these_mice=[these_mice handles.mouse_no(fileNo)];
% 
% end
% 
% pffft=1;
% 
% end
% 
% these_delta_log_P_laser=these_delta_log_P_laser(~isnan(these_delta_log_P_laser));
% 
% these_mice_laser=these_mice(~isnan(these_delta_log_P_laser));
% 
% these_delta_log_P_post=these_delta_log_P_post(~isnan(these_delta_log_P_post));
% 
% these_mice_post=these_mice(~isnan(these_delta_log_P_post));
% 
% %Laser
% 
% bar(bar_offset,mean(these_delta_log_P_laser),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
% 
% [mean_out, CIout]=drgViolinPoint(these_delta_log_P_laser,edges,bar_offset,rand_offset,'k','k',3);
% 
% glm_lfp_power.data(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=these_delta_log_P_laser;
% 
% glm_lfp_power.time_window(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=0*ones(1,length(these_delta_log_P_laser));
% 
% switch groupNo
% 
% case 1
% 
% glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=1;
% 
% glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=1;
% 
% case 2
% 
% glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=1;
% 
% glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=2;
% 
% case 3
% 
% glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=2;
% 
% glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=1;
% 
% case 4
% 
% glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=2;
% 
% glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=2;
% 
% end
% 
% glm_lfp_power.mouseNo(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=these_mice_laser;
% 
% glm_ii_lfp=glm_ii_lfp+length(these_delta_log_P_laser);
% 
% id_ii=id_ii+1;
% 
% input_data(id_ii).data=these_delta_log_P_laser;
% 
% input_data(id_ii).description=['laser ' group_label{groupNo}];
% 
% %Post
% 
% bar_offset=bar_offset+1;
% 
% bar(bar_offset,mean(these_delta_log_P_post),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
% 
% [mean_out, CIout]=drgViolinPoint(these_delta_log_P_post,edges,bar_offset,rand_offset,'k','k',3);
% 
% glm_lfp_power.data(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=these_delta_log_P_post;
% 
% glm_lfp_power.time_window(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=1*ones(1,length(these_delta_log_P_post));
% 
% switch groupNo
% 
% case 1
% 
% glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=1;
% 
% glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=1;
% 
% case 2
% 
% glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=1;
% 
% glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=2;
% 
% case 3
% 
% glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=2;
% 
% glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=1;
% 
% case 4
% 
% glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=2;
% 
% glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=2;
% 
% end
% 
% glm_lfp_power.mouseNo(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=these_mice_post;
% 
% glm_ii_lfp=glm_ii_lfp+length(these_delta_log_P_post);
% 
% id_ii=id_ii+1;
% 
% input_data(id_ii).data=these_delta_log_P_post;
% 
% input_data(id_ii).description=['post ' group_label{groupNo}];
% 
% bar_offset=bar_offset+2;
% 
% else
% 
% bar_offset=bar_offset+3;
% 
% end
% 
% end
% 
% title(['Peak Angle Variance ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}])
% 
% xticks([0.5 3.5 6.5 9.5])
% 
% xticklabels({group_label{1},group_label{2},group_label{3},group_label{4}})
% 
% % ylim([0 0.03])
% 
% ylabel('Angle Variance')
% 
% %Perform the glm for MI per mouse per odor pair
% 
% fprintf(1, ['glm for peak angle variance for ' handles.electrode_label{handles.peakLFPNo(ii_electrode)} '\n'])
% 
% fprintf(fileID, ['glm for peak angle variance for ' handles.electrode_label{handles.peakLFPNo(ii_electrode)} '\n']);
% 
% % tbl = table(glm_lfp_power.data',glm_lfp_power.time_window',glm_lfp_power.group',...
% 
% % 'VariableNames',{'logP','laser_post','group'});
% 
% % mdl = fitglm(tbl,'logP~laser_post+group+laser_post*group'...
% 
% % ,'CategoricalVars',[2,3])
% 
% % First, create a cell array of character labels corresponding to your group numbers
% 
% % group_labels = {'Control', '5xFAD', 'Group3', 'Group4'}; % Add all your group labels here
% 
% tbl = table(glm_lfp_power.data',categorical(glm_lfp_power.time_window', [0 1], time_window_labels),categorical(glm_lfp_power.genotype',[1], genotype_label),categorical(glm_lfp_power.laser',[1], laser_label),'VariableNames',{'logP','time_window','genotype','laser'});
% 
% mdl = fitglm(tbl,'logP~time_window+genotype+laser+time_window*genotype*laser');
% 
% txt = evalc('mdl');
% 
% txt=regexp(txt,'','split');
% 
% txt=cell2mat(txt);
% 
% txt=regexp(txt,'','split');
% 
% txt=cell2mat(txt);
% 
% fprintf(fileID,'%s\n', txt);
% 
% %Do the ranksum/t-test
% 
% fprintf(1, ['\n\nRanksum or t-test p values for peak angle variance for ' handles.electrode_label{handles.peakLFPNo(ii_electrode)} '\n']);
% 
% fprintf(fileID, ['\n\nRanksum or t-test p values for peak angle variance for ' handles.electrode_label{handles.peakLFPNo(ii_electrode)} '\n']);
% 
% [output_data] = drgMutiRanksumorTtest(input_data,fileID);



end
% toc
% fclose(fileID);



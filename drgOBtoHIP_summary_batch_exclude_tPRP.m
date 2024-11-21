function drgOBtoHIP_summary_batch_exclude_tPRP(choiceBatchPathName,choiceFileName)
%Note: fitcnet will not work in Matlab versions earlier than 2021a

close all
clear all

PathName_out='/Users/restrepd/Documents/Projects/Joe_OB_to_hippo/5xFADvsWT_1_hour_treatmentVsnone/OBtoHIP_out/';

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


our_colors(1).color=[1 0 0];
our_colors(2).color=[0 1 0];
our_colors(3).color=[0 0 1];
our_colors(4).color=[1 1 0];

%Re-arrange the groups
%Group numbers
%1 WT Treated
%2 5xFAD Treated
%3 WT Untreated
%4 5xFAD Untreated


group_label{1}='WT U';
group_label{2}='WT T';
group_label{3}='5xFAD U';
group_label{4}='5xFAD T';

genotype_label{1}='WT';
genotype_label{2}='5xFAD';

laser_label{1}='UT';
laser_label{2}='T';

time_window_labels = {'Laser', 'Post'}; % Add all your group labels here
logP_time_window_labels = {'Pre','Laser', 'Post'};

t_ref_start=0;
t_ref_end=5*60; %In seconds
t_laser_end=65*60;


%You can choose which plots are shown
these_plots_shown(1)=1; %delta peak tPRP bar graph for laser and post windows (minus pre)
these_plots_shown(2)=1; %delta peak tPRP bar graph for laser window (minus pre)
these_plots_shown(3)=0; %delta trough tPRP bar graph for laser and post windows (minus pre)
these_plots_shown(4)=1; %delta trough tPRP bar graph for laser window (minus pre)
these_plots_shown(5)=0; %tPRP peak for each electrode and each time window
these_plots_shown(6)=1; %tPRP peak for each electrode for laser window only
these_plots_shown(7)=0; %tPRP trough for each electrode and each time window
these_plots_shown(8)=1; %tPRP trough for each electrode for laser window only
these_plots_shown(9)=0; %show modulation index for each electrode for each time window
these_plots_shown(10)=1; %show modulation index for each electrode for each time window
these_plots_shown(11)=0; %show modulation index as a function of session number for each electrode for only one time window
these_plots_shown(12)=0; %show peak angle variance for each electrode for each time window
these_plots_shown(13)=1; %show peak angle variance as a function of session number for each electrode for only one time window

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

    %Show delta peak tPRP bar graphs for each electrode for both laser and post
    if these_plots_shown(1)==1
        edges=[-30:5:15];
        rand_offset=0.8;


        for ii_electrode=1:4
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            hFig2 = figure(figNo);
            set(hFig2, 'units','normalized','position',[.1 .1 .3 .3])
            hold on

            bar_offset=0;

            glm_lfp_power=[];
            glm_ii_lfp=0;

            id_ii=0;
            input_data=[];

            for groupNo=1:length(groups)
                %Find the files for this mouse
                these_files=[];

                for fileNo=1:handles.no_files
                    if handles.group_no(fileNo)==groupNo
                        %Does this file exist?
                        this_filename=handles.FileName{fileNo};
                        FileName_out=['OBtoHIPtPRP_' this_filename(11:end)];
                        if exist([PathName_out FileName_out])~=0
                            these_files=[these_files fileNo];
                        end
                    end
                end

                %Get deltaP
                these_delta_log_P_laser=[];
                these_delta_log_P_post=[];
                these_mice=[];

                if ~isempty(these_files)
                    for fileNo=these_files
                        %Plot the timecourse
                        % PathName_out=handles.PathName{fileNo};
                        this_filename=handles.FileName{fileNo};
                        FileName_out=['OBtoHIPtPRP_' this_filename(11:end)];
                        load([PathName_out FileName_out])
                        add_electrode=0;
                        if isempty(handles.file(fileNo).exclude_electrodes)
                            add_electrode=1;
                        else
                            if sum(handles.file(fileNo).exclude_electrodes==actual_electrode(ii_electrode))==0
                                add_electrode=1;
                            end

                        end
                        if add_electrode==1
                            log_P_timecourse=handles_out.electrode(ii_electrode).handles_out.mean_peakPower;
                            time=handles_out.electrode(ii_electrode).handles_out.time;
                            this_delta_log_P_laser=mean(mean(log_P_timecourse((time>=t_ref_end)&(time<=t_laser_end))))-...
                                mean(mean(log_P_timecourse(time<=t_ref_end)));
                            these_delta_log_P_laser=[these_delta_log_P_laser this_delta_log_P_laser];
                            this_delta_log_P_post=mean(mean(log_P_timecourse(time>=t_laser_end)))-...
                                mean(mean(log_P_timecourse(time<=t_ref_end)));
                            these_delta_log_P_post=[these_delta_log_P_post this_delta_log_P_post];
                            these_mice=[these_mice handles.mouse_no(fileNo)];
                        end
                        pffft=1;
                    end

                    these_delta_log_P_laser=these_delta_log_P_laser(~isnan(these_delta_log_P_laser));
                    these_mice_laser=these_mice(~isnan(these_delta_log_P_laser));
                    these_delta_log_P_post=these_delta_log_P_post(~isnan(these_delta_log_P_post));
                    these_mice_post=these_mice(~isnan(these_delta_log_P_post));


                    %Laser
                    bar(bar_offset,mean(these_delta_log_P_laser),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
                    [mean_out, CIout]=drgViolinPoint(these_delta_log_P_laser,edges,bar_offset,rand_offset,'k','k',3);

                    glm_lfp_power.data(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=these_delta_log_P_laser;
                    glm_lfp_power.time_window(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=0*ones(1,length(these_delta_log_P_laser));
                    switch groupNo
                        case 1
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=1;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=1;
                        case 2
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=1;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=2;
                        case 3
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=2;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=1;
                        case 4
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=2;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=2;
                    end
                    glm_lfp_power.mouseNo(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=these_mice_laser;
                    glm_ii_lfp=glm_ii_lfp+length(these_delta_log_P_laser);

                    id_ii=id_ii+1;
                    input_data(id_ii).data=these_delta_log_P_laser;
                    input_data(id_ii).description=['laser ' group_label{groupNo}];

                    %Post
                    bar_offset=bar_offset+1;
                    bar(bar_offset,mean(these_delta_log_P_post),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])

                    [mean_out, CIout]=drgViolinPoint(these_delta_log_P_post,edges,bar_offset,rand_offset,'k','k',3);

                    glm_lfp_power.data(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=these_delta_log_P_post;
                    glm_lfp_power.time_window(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=1*ones(1,length(these_delta_log_P_post));
                    switch groupNo
                        case 1
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=1;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=1;
                        case 2
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=1;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=2;
                        case 3
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=2;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=1;
                        case 4
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=2;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=2;
                    end
                    glm_lfp_power.mouseNo(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=these_mice_post;
                    glm_ii_lfp=glm_ii_lfp+length(these_delta_log_P_post);

                    id_ii=id_ii+1;
                    input_data(id_ii).data=these_delta_log_P_post;
                    input_data(id_ii).description=['post ' group_label{groupNo}];

                    bar_offset=bar_offset+2;
                else
                    bar_offset=bar_offset+3;
                end
            end

            title(['delta peak tPRP ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}])


            xticks([0.5 3.5 6.5 9.5])
            xticklabels({group_label{1}, group_label{2}....
                ,group_label{3}, group_label{4}})


            % ylim([0 0.03])

            ylabel('dB')

            %Perform the glm for MI per mouse per odor pair
            fprintf(1, ['glm for delta peak tPRP for '  handles.electrode_label{handles.peakLFPNo(ii_electrode)} '\n'])
            fprintf(fileID, ['glm for delta peak tPRP for '  handles.electrode_label{handles.peakLFPNo(ii_electrode)} '\n']);

            % tbl = table(glm_lfp_power.data',glm_lfp_power.time_window',glm_lfp_power.group',...
            %     'VariableNames',{'logP','laser_post','group'});
            % mdl = fitglm(tbl,'logP~laser_post+group+laser_post*group'...
            %     ,'CategoricalVars',[2,3])

            % First, create a cell array of character labels corresponding to your group numbers
            % group_labels = {'Control', '5xFAD', 'Group3', 'Group4'}; % Add all your group labels here

            tbl = table(glm_lfp_power.data',...
                categorical(glm_lfp_power.time_window', [0 1], time_window_labels),...
                categorical(glm_lfp_power.genotype', [1 2], genotype_label),...
                categorical(glm_lfp_power.laser', [1 2], laser_label),...
                'VariableNames',{'logP','time_window','genotype','laser'});
            mdl = fitglm(tbl,'logP~time_window+genotype+laser+time_window*genotype*laser')

            txt = evalc('mdl');
            txt=regexp(txt,'<strong>','split');
            txt=cell2mat(txt);
            txt=regexp(txt,'</strong>','split');
            txt=cell2mat(txt);

            fprintf(fileID,'%s\n', txt);


            %Do the ranksum/t-test
            fprintf(1, ['\n\nRanksum or t-test p values for delta peak tPRP for ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}  '\n']);
            fprintf(fileID, ['\n\nRanksum or t-test p values for delta peak tPRP for ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}  '\n']);

            [output_data] = drgMutiRanksumorTtest(input_data,fileID);

        end
    end


    %Show delta peak tPRP bar graphs for each electrode for laser
    if these_plots_shown(2)==1
        edges=[-30:5:15];
        rand_offset=0.8;


        for ii_electrode=1:4
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            hFig2 = figure(figNo);
            set(hFig2, 'units','normalized','position',[.1 .1 .3 .3])
            hold on

            bar_offset=0;

            glm_lfp_power=[];
            glm_ii_lfp=0;

            id_ii=0;
            input_data=[];

            for groupNo=1:length(groups)
                %Find the files for this mouse
                these_files=[];

                for fileNo=1:handles.no_files
                    if handles.group_no(fileNo)==groupNo
                        %Does this file exist?
                        this_filename=handles.FileName{fileNo};
                        FileName_out=['OBtoHIPtPRP_' this_filename(11:end)];
                        if exist([PathName_out FileName_out])~=0
                            these_files=[these_files fileNo];
                        end
                    end
                end

                %Get deltaP
                these_delta_log_P_laser=[];
                these_delta_log_P_post=[];
                these_mice=[];

                if ~isempty(these_files)
                    for fileNo=these_files
                        %Plot the timecourse
                        % PathName_out=handles.PathName{fileNo};
                        this_filename=handles.FileName{fileNo};
                        FileName_out=['OBtoHIPtPRP_' this_filename(11:end)];
                        load([PathName_out FileName_out])
                        add_electrode=0;
                        if isempty(handles.file(fileNo).exclude_electrodes)
                            add_electrode=1;
                        else
                            if sum(handles.file(fileNo).exclude_electrodes==actual_electrode(ii_electrode))==0
                                add_electrode=1;
                            end

                        end
                        if add_electrode==1
                            log_P_timecourse=handles_out.electrode(ii_electrode).handles_out.mean_peakPower;
                            time=handles_out.electrode(ii_electrode).handles_out.time;
                            this_delta_log_P_laser=mean(mean(log_P_timecourse((time>=t_ref_end)&(time<=t_laser_end))))-...
                                mean(mean(log_P_timecourse(time<=t_ref_end)));
                            these_delta_log_P_laser=[these_delta_log_P_laser this_delta_log_P_laser];
                            this_delta_log_P_post=mean(mean(log_P_timecourse(time>=t_laser_end)))-...
                                mean(mean(log_P_timecourse(time<=t_ref_end)));
                            these_delta_log_P_post=[these_delta_log_P_post this_delta_log_P_post];
                            these_mice=[these_mice handles.mouse_no(fileNo)];
                        end
                        pffft=1;
                    end

                    these_delta_log_P_laser=these_delta_log_P_laser(~isnan(these_delta_log_P_laser));
                    these_mice_laser=these_mice(~isnan(these_delta_log_P_laser));
                    these_delta_log_P_post=these_delta_log_P_post(~isnan(these_delta_log_P_post));
                    these_mice_post=these_mice(~isnan(these_delta_log_P_post));


                    %Laser
                    bar(bar_offset,mean(these_delta_log_P_laser),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
                    [mean_out, CIout]=drgViolinPoint(these_delta_log_P_laser,edges,bar_offset,rand_offset,'k','k',3);

                    glm_lfp_power.data(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=these_delta_log_P_laser;
                    glm_lfp_power.time_window(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=0*ones(1,length(these_delta_log_P_laser));
                    switch groupNo
                        case 1
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=1;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=1;
                        case 2
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=1;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=2;
                        case 3
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=2;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=1;
                        case 4
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=2;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=2;
                    end
                    glm_lfp_power.mouseNo(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=these_mice_laser;
                    glm_ii_lfp=glm_ii_lfp+length(these_delta_log_P_laser);

                    id_ii=id_ii+1;
                    input_data(id_ii).data=these_delta_log_P_laser;
                    input_data(id_ii).description=['laser ' group_label{groupNo}];

                    bar_offset=bar_offset+1;
                else
                    bar_offset=bar_offset+3;
                end
            end

            title(['delta peak tPRP ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}])


            xticks([0 1 2 3])
            xticklabels({group_label{1}, group_label{2}....
                ,group_label{3}, group_label{4}})


            % ylim([0 0.03])

            ylabel('dB')

            %Perform the glm for MI per mouse per odor pair
            fprintf(1, ['glm for delta peak tPRP for '  handles.electrode_label{handles.peakLFPNo(ii_electrode)} '\n'])
            fprintf(fileID, ['glm for delta peak tPRP for '  handles.electrode_label{handles.peakLFPNo(ii_electrode)} '\n']);

            % tbl = table(glm_lfp_power.data',glm_lfp_power.time_window',glm_lfp_power.group',...
            %     'VariableNames',{'logP','laser_post','group'});
            % mdl = fitglm(tbl,'logP~laser_post+group+laser_post*group'...
            %     ,'CategoricalVars',[2,3])

            % First, create a cell array of character labels corresponding to your group numbers
            % group_labels = {'Control', '5xFAD', 'Group3', 'Group4'}; % Add all your group labels here

            tbl = table(glm_lfp_power.data',...
                categorical(glm_lfp_power.genotype', [1 2], genotype_label),...
                categorical(glm_lfp_power.laser', [1 2], laser_label),...
                'VariableNames',{'logP','genotype','laser'});
            mdl = fitglm(tbl,'logP~genotype+laser+genotype*laser')

            txt = evalc('mdl');
            txt=regexp(txt,'<strong>','split');
            txt=cell2mat(txt);
            txt=regexp(txt,'</strong>','split');
            txt=cell2mat(txt);

            fprintf(fileID,'%s\n', txt);


            %Do the ranksum/t-test
            fprintf(1, ['\n\nRanksum or t-test p values for delta peak tPRP for ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}  '\n']);
            fprintf(fileID, ['\n\nRanksum or t-test p values for delta peak tPRP for ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}  '\n']);

            [output_data] = drgMutiRanksumorTtest(input_data,fileID);

        end
    end
    

    %Show delta trough tPRP bar graphs for each electrode for both laser and
    %post
    if these_plots_shown(3)==1
        for ii_electrode=1:4
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            hFig2 = figure(figNo);
            set(hFig2, 'units','normalized','position',[.1 .1 .3 .3])
            hold on

            bar_offset=0;

            glm_lfp_power=[];
            glm_ii_lfp=0;

            id_ii=0;
            input_data=[];

            for groupNo=1:length(groups)
                %Find the files for this mouse
                these_files=[];

                for fileNo=1:handles.no_files
                    if handles.group_no(fileNo)==groupNo
                        %Does this file exist?
                        this_filename=handles.FileName{fileNo};
                        FileName_out=['OBtoHIPtPRP_' this_filename(11:end)];
                        if exist([PathName_out FileName_out])~=0
                            these_files=[these_files fileNo];
                        end
                    end
                end

                %Get deltaP
                these_delta_log_P_laser=[];
                these_delta_log_P_post=[];
                these_mice=[];

                if ~isempty(these_files)
                    for fileNo=these_files
                        %Plot the timecourse
                        % PathName_out=handles.PathName{fileNo};
                        this_filename=handles.FileName{fileNo};
                        FileName_out=['OBtoHIPtPRP_' this_filename(11:end)];
                        load([PathName_out FileName_out])
                        add_electrode=0;
                        if isempty(handles.file(fileNo).exclude_electrodes)
                            add_electrode=1;
                        else
                            if sum(handles.file(fileNo).exclude_electrodes==actual_electrode(ii_electrode))==0
                                add_electrode=1;
                            end

                        end
                        if add_electrode==1
                            log_P_timecourse=handles_out.electrode(ii_electrode).handles_out.mean_troughPower;
                            time=handles_out.electrode(ii_electrode).handles_out.time;
                            this_delta_log_P_laser=mean(mean(log_P_timecourse((time>=t_ref_end)&(time<=t_laser_end))))-...
                                mean(mean(log_P_timecourse(time<=t_ref_end)));
                            these_delta_log_P_laser=[these_delta_log_P_laser this_delta_log_P_laser];
                            this_delta_log_P_post=mean(mean(log_P_timecourse(time>=t_laser_end)))-...
                                mean(mean(log_P_timecourse(time<=t_ref_end)));
                            these_delta_log_P_post=[these_delta_log_P_post this_delta_log_P_post];
                            these_mice=[these_mice handles.mouse_no(fileNo)];
                        end
                        pffft=1;
                    end

                    these_delta_log_P_laser=these_delta_log_P_laser(~isnan(these_delta_log_P_laser));
                    these_mice_laser=these_mice(~isnan(these_delta_log_P_laser));
                    these_delta_log_P_post=these_delta_log_P_post(~isnan(these_delta_log_P_post));
                    these_mice_post=these_mice(~isnan(these_delta_log_P_post));


                    %Laser
                    bar(bar_offset,mean(these_delta_log_P_laser),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
                    [mean_out, CIout]=drgViolinPoint(these_delta_log_P_laser,edges,bar_offset,rand_offset,'k','k',3);

                    glm_lfp_power.data(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=these_delta_log_P_laser;
                    glm_lfp_power.time_window(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=0*ones(1,length(these_delta_log_P_laser));
                    switch groupNo
                        case 1
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=1;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=1;
                        case 2
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=1;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=2;
                        case 3
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=2;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=1;
                        case 4
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=2;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=2;
                    end
                    glm_lfp_power.mouseNo(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=these_mice_laser;
                    glm_ii_lfp=glm_ii_lfp+length(these_delta_log_P_laser);

                    id_ii=id_ii+1;
                    input_data(id_ii).data=these_delta_log_P_laser;
                    input_data(id_ii).description=['laser ' group_label{groupNo}];

                    %Post
                    bar_offset=bar_offset+1;
                    bar(bar_offset,mean(these_delta_log_P_post),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])

                    [mean_out, CIout]=drgViolinPoint(these_delta_log_P_post,edges,bar_offset,rand_offset,'k','k',3);

                    glm_lfp_power.data(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=these_delta_log_P_post;
                    glm_lfp_power.time_window(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=1*ones(1,length(these_delta_log_P_post));
                    switch groupNo
                        case 1
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=1;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=1;
                        case 2
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=1;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=2;
                        case 3
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=2;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=1;
                        case 4
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=2;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=2;
                    end
                    glm_lfp_power.mouseNo(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=these_mice_post;
                    glm_ii_lfp=glm_ii_lfp+length(these_delta_log_P_post);

                    id_ii=id_ii+1;
                    input_data(id_ii).data=these_delta_log_P_post;
                    input_data(id_ii).description=['post ' group_label{groupNo}];

                    bar_offset=bar_offset+2;
                else
                    bar_offset=bar_offset+3;
                end
            end

            title(['delta trough tPRP ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}])


            xticks([0.5 3.5 6.5 9.5])
            xticklabels({group_label{1}, group_label{2}....
                ,group_label{3}, group_label{4}})


            % ylim([0 0.03])

            ylabel('dB')

            %Perform the glm for MI per mouse per odor pair
            fprintf(1, ['glm for delta trough tPRP for '  handles.electrode_label{handles.peakLFPNo(ii_electrode)} '\n'])
            fprintf(fileID, ['glm for delta trough tPRP for '  handles.electrode_label{handles.peakLFPNo(ii_electrode)} '\n']);

            % tbl = table(glm_lfp_power.data',glm_lfp_power.time_window',glm_lfp_power.group',...
            %     'VariableNames',{'logP','laser_post','group'});
            % mdl = fitglm(tbl,'logP~laser_post+group+laser_post*group'...
            %     ,'CategoricalVars',[2,3])

            % First, create a cell array of character labels corresponding to your group numbers
            % group_labels = {'Control', '5xFAD', 'Group3', 'Group4'}; % Add all your group labels here

            tbl = table(glm_lfp_power.data',...
                categorical(glm_lfp_power.time_window', [0 1], time_window_labels),...
                categorical(glm_lfp_power.genotype', [1 2], genotype_label),...
                categorical(glm_lfp_power.laser', [1 2], laser_label),...
                'VariableNames',{'logP','time_window','genotype','laser'});
            mdl = fitglm(tbl,'logP~time_window+genotype+laser+time_window*genotype*laser')

            txt = evalc('mdl');
            txt=regexp(txt,'<strong>','split');
            txt=cell2mat(txt);
            txt=regexp(txt,'</strong>','split');
            txt=cell2mat(txt);

            fprintf(fileID,'%s\n', txt);


            %Do the ranksum/t-test
            fprintf(1, ['\n\nRanksum or t-test p values for delta trough tPRP for ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}  '\n']);
            fprintf(fileID, ['\n\nRanksum or t-test p values for delta trough tPRP for ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}  '\n']);

            [output_data] = drgMutiRanksumorTtest(input_data,fileID);

        end
    end

    %Show delta trough tPRP bar graphs for each electrode for laser only
    if these_plots_shown(4)==1
        for ii_electrode=1:4
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            hFig2 = figure(figNo);
            set(hFig2, 'units','normalized','position',[.1 .1 .3 .3])
            hold on

            bar_offset=0;

            glm_lfp_power=[];
            glm_ii_lfp=0;

            id_ii=0;
            input_data=[];

            for groupNo=1:length(groups)
                %Find the files for this mouse
                these_files=[];

                for fileNo=1:handles.no_files
                    if handles.group_no(fileNo)==groupNo
                        %Does this file exist?
                        this_filename=handles.FileName{fileNo};
                        FileName_out=['OBtoHIPtPRP_' this_filename(11:end)];
                        if exist([PathName_out FileName_out])~=0
                            these_files=[these_files fileNo];
                        end
                    end
                end

                %Get deltaP
                these_delta_log_P_laser=[];
                these_delta_log_P_post=[];
                these_mice=[];

                if ~isempty(these_files)
                    for fileNo=these_files
                        %Plot the timecourse
                        % PathName_out=handles.PathName{fileNo};
                        this_filename=handles.FileName{fileNo};
                        FileName_out=['OBtoHIPtPRP_' this_filename(11:end)];
                        load([PathName_out FileName_out])
                        add_electrode=0;
                        if isempty(handles.file(fileNo).exclude_electrodes)
                            add_electrode=1;
                        else
                            if sum(handles.file(fileNo).exclude_electrodes==actual_electrode(ii_electrode))==0
                                add_electrode=1;
                            end

                        end
                        if add_electrode==1
                            log_P_timecourse=handles_out.electrode(ii_electrode).handles_out.mean_troughPower;
                            time=handles_out.electrode(ii_electrode).handles_out.time;
                            this_delta_log_P_laser=mean(mean(log_P_timecourse((time>=t_ref_end)&(time<=t_laser_end))))-...
                                mean(mean(log_P_timecourse(time<=t_ref_end)));
                            these_delta_log_P_laser=[these_delta_log_P_laser this_delta_log_P_laser];
                            this_delta_log_P_post=mean(mean(log_P_timecourse(time>=t_laser_end)))-...
                                mean(mean(log_P_timecourse(time<=t_ref_end)));
                            these_delta_log_P_post=[these_delta_log_P_post this_delta_log_P_post];
                            these_mice=[these_mice handles.mouse_no(fileNo)];
                        end
                        pffft=1;
                    end

                    these_delta_log_P_laser=these_delta_log_P_laser(~isnan(these_delta_log_P_laser));
                    these_mice_laser=these_mice(~isnan(these_delta_log_P_laser));
                    these_delta_log_P_post=these_delta_log_P_post(~isnan(these_delta_log_P_post));
                    these_mice_post=these_mice(~isnan(these_delta_log_P_post));


                    %Laser
                    bar(bar_offset,mean(these_delta_log_P_laser),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
                    [mean_out, CIout]=drgViolinPoint(these_delta_log_P_laser,edges,bar_offset,rand_offset,'k','k',3);

                    glm_lfp_power.data(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=these_delta_log_P_laser;
                    glm_lfp_power.time_window(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=0*ones(1,length(these_delta_log_P_laser));
                    switch groupNo
                        case 1
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=1;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=1;
                        case 2
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=1;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=2;
                        case 3
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=2;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=1;
                        case 4
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=2;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=2;
                    end
                    glm_lfp_power.mouseNo(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=these_mice_laser;
                    glm_ii_lfp=glm_ii_lfp+length(these_delta_log_P_laser);

                    id_ii=id_ii+1;
                    input_data(id_ii).data=these_delta_log_P_laser;
                    input_data(id_ii).description=['laser ' group_label{groupNo}];


                    bar_offset=bar_offset+1;
                else
                    bar_offset=bar_offset+3;
                end
            end

            title(['delta trough tPRP for laser ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}])


            xticks([0 1 2 3])
            xticklabels({group_label{1}, group_label{2}....
                ,group_label{3}, group_label{4}})


            % ylim([0 0.03])

            ylabel('dB')

            %Perform the glm for MI per mouse per odor pair
            fprintf(1, ['glm for delta trough tPRP for laser for '  handles.electrode_label{handles.peakLFPNo(ii_electrode)} '\n'])
            fprintf(fileID, ['glm for delta trough tPRP for laser for '  handles.electrode_label{handles.peakLFPNo(ii_electrode)} '\n']);

            % tbl = table(glm_lfp_power.data',glm_lfp_power.time_window',glm_lfp_power.group',...
            %     'VariableNames',{'logP','laser_post','group'});
            % mdl = fitglm(tbl,'logP~laser_post+group+laser_post*group'...
            %     ,'CategoricalVars',[2,3])

            % First, create a cell array of character labels corresponding to your group numbers
            % group_labels = {'Control', '5xFAD', 'Group3', 'Group4'}; % Add all your group labels here

            tbl = table(glm_lfp_power.data',...
                categorical(glm_lfp_power.genotype', [1 2], genotype_label),...
                categorical(glm_lfp_power.laser', [1 2], laser_label),...
                'VariableNames',{'logP','genotype','laser'});
            mdl = fitglm(tbl,'logP~genotype+laser+genotype*laser')

            txt = evalc('mdl');
            txt=regexp(txt,'<strong>','split');
            txt=cell2mat(txt);
            txt=regexp(txt,'</strong>','split');
            txt=cell2mat(txt);

            fprintf(fileID,'%s\n', txt);


            %Do the ranksum/t-test
            fprintf(1, ['\n\nRanksum or t-test p values for delta trough tPRP for laser ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}  '\n']);
            fprintf(fileID, ['\n\nRanksum or t-test p values for delta trough tPRP for laser ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}  '\n']);

            [output_data] = drgMutiRanksumorTtest(input_data,fileID);

        end
    end






    edges=[10:5:80];
    rand_offset=0.8;

    %Now show tPRP peak for each electrode and each time window
    if these_plots_shown(5)==1
        for ii_electrode=1:4
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            hFig2 = figure(figNo);
            set(hFig2, 'units','normalized','position',[.1 .1 .3 .3])
            hold on

            bar_offset=0;

            glm_lfp_power=[];
            glm_ii_lfp=0;

            id_ii=0;
            input_data=[];

            for groupNo=1:length(groups)
                %Find the files for this mouse
                these_files=[];

                for fileNo=1:handles.no_files
                    if handles.group_no(fileNo)==groupNo
                        %Does this file exist?
                        this_filename=handles.FileName{fileNo};
                        FileName_out=['OBtoHIPtPRP_' this_filename(11:end)];
                        if exist([PathName_out FileName_out])~=0
                            these_files=[these_files fileNo];
                        end
                    end
                end

                %Get deltaP
                these_log_P_pre=[];
                these_log_P_laser=[];
                these_log_P_post=[];
                these_mice=[];

                if ~isempty(these_files)
                    for fileNo=these_files
                        %Plot the timecourse
                        % PathName_out=handles.PathName{fileNo};
                        this_filename=handles.FileName{fileNo};
                        FileName_out=['OBtoHIPtPRP_' this_filename(11:end)];
                        load([PathName_out FileName_out])
                        add_electrode=0;
                        if isempty(handles.file(fileNo).exclude_electrodes)
                            add_electrode=1;
                        else
                            if sum(handles.file(fileNo).exclude_electrodes==actual_electrode(ii_electrode))==0
                                add_electrode=1;
                            end

                        end
                        if add_electrode==1
                            log_P_timecourse=handles_out.electrode(ii_electrode).handles_out.mean_peakPower;
                            time=handles_out.electrode(ii_electrode).handles_out.time;
                            this_log_P_pre=mean(mean(log_P_timecourse(time<t_ref_end)));
                            these_log_P_pre=[these_log_P_pre this_log_P_pre];
                            this_log_P_laser=mean(mean(log_P_timecourse((time>=t_ref_end)&(time<=t_laser_end))));
                            these_log_P_laser=[these_log_P_laser this_log_P_laser];
                            this_log_P_post=mean(mean(log_P_timecourse(time>t_laser_end)));
                            these_log_P_post=[these_log_P_post this_log_P_post];
                            these_mice=[these_mice handles.mouse_no(fileNo)];
                        end
                        pffft=1;
                    end

                    these_log_P_pre=these_log_P_pre(~isnan(these_log_P_pre));
                    these_mice_pre=these_mice(~isnan(these_log_P_pre));
                    these_log_P_laser=these_log_P_laser(~isnan(these_log_P_laser));
                    these_mice_laser=these_mice(~isnan(these_log_P_laser));
                    these_log_P_post=these_log_P_post(~isnan(these_log_P_post));
                    these_mice_post=these_mice(~isnan(these_log_P_post));


                    %Pre
                    bar(bar_offset,mean(these_log_P_pre),'LineWidth', 3,'EdgeColor','none','FaceColor',[31/255 31/255 140/255])
                    [mean_out, CIout]=drgViolinPoint(these_log_P_pre,edges,bar_offset,rand_offset,'k','k',3);

                    glm_lfp_power.data(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=these_log_P_pre;
                    glm_lfp_power.time_window(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=0*ones(1,length(these_log_P_pre));
                    switch groupNo
                        case 1
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=1;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=1;
                        case 2
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=1;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=2;
                        case 3
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=2;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=1;
                        case 4
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=2;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=2;
                    end
                    glm_lfp_power.mouseNo(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=these_mice_pre;
                    glm_ii_lfp=glm_ii_lfp+length(these_log_P_pre);

                    id_ii=id_ii+1;
                    input_data(id_ii).data=these_log_P_pre;
                    input_data(id_ii).description=['laser ' group_label{groupNo}];

                    %Laser
                    bar_offset=bar_offset+1;
                    bar(bar_offset,mean(these_log_P_laser),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
                    [mean_out, CIout]=drgViolinPoint(these_log_P_laser,edges,bar_offset,rand_offset,'k','k',3);

                    glm_lfp_power.data(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=these_log_P_laser;
                    glm_lfp_power.time_window(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=1*ones(1,length(these_log_P_laser));
                    switch groupNo
                        case 1
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=1;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=1;
                        case 2
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=1;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=2;
                        case 3
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=2;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=1;
                        case 4
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=2;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=2;
                    end
                    glm_lfp_power.mouseNo(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=these_mice_laser;
                    glm_ii_lfp=glm_ii_lfp+length(these_log_P_laser);

                    id_ii=id_ii+1;
                    input_data(id_ii).data=these_log_P_laser;
                    input_data(id_ii).description=['laser ' group_label{groupNo}];

                    %Post
                    bar_offset=bar_offset+1;
                    bar(bar_offset,mean(these_log_P_post),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
                    [mean_out, CIout]=drgViolinPoint(these_log_P_post,edges,bar_offset,rand_offset,'k','k',3);

                    glm_lfp_power.data(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=these_log_P_post;
                    glm_lfp_power.time_window(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=2*ones(1,length(these_log_P_post));
                    switch groupNo
                        case 1
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=1;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=1;
                        case 2
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=1;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=2;
                        case 3
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=2;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=1;
                        case 4
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=2;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=2;
                    end
                    glm_lfp_power.mouseNo(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=these_mice_post;
                    glm_ii_lfp=glm_ii_lfp+length(these_log_P_post);

                    id_ii=id_ii+1;
                    input_data(id_ii).data=these_log_P_post;
                    input_data(id_ii).description=['post ' group_label{groupNo}];

                    bar_offset=bar_offset+2;
                else
                    bar_offset=bar_offset+3;
                end
            end

            title(['Peak tPRP ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}])


            xticks([1 5 9 13])
            xticklabels({group_label{1}, group_label{2}....
                ,group_label{3}, group_label{4}})


            % ylim([0 0.03])

            ylabel('dB')

            %Perform the glm for MI per mouse per odor pair
            fprintf(1, ['glm for peak tPRP for ' handles.electrode_label{handles.peakLFPNo(ii_electrode)} '\n'])
            fprintf(fileID, ['glm for peak tPRP for '  handles.electrode_label{handles.peakLFPNo(ii_electrode)} '\n']);

            % tbl = table(glm_lfp_power.data',glm_lfp_power.time_window',glm_lfp_power.group',...
            %     'VariableNames',{'logP','laser_post','group'});
            % mdl = fitglm(tbl,'logP~laser_post+group+laser_post*group'...
            %     ,'CategoricalVars',[2,3])

            % First, create a cell array of character labels corresponding to your group numbers
            % group_labels = {'Control', '5xFAD', 'Group3', 'Group4'}; % Add all your group labels here


            % for ii=1:length(re_arranged_groups)
            %     re_arranged_group_labels{ii}=handles.group_label{re_arranged_groups(ii)};
            % end
            tbl = table(glm_lfp_power.data',...
                categorical(glm_lfp_power.time_window', [0 1 2], logP_time_window_labels),...
                categorical(glm_lfp_power.genotype', [1 2], genotype_label),...
                categorical(glm_lfp_power.laser', [1 2], laser_label),...
                'VariableNames',{'logP','time_window','genotype','laser'});
            mdl = fitglm(tbl,'logP~time_window+genotype+laser+time_window*genotype*laser')

            txt = evalc('mdl');
            txt=regexp(txt,'<strong>','split');
            txt=cell2mat(txt);
            txt=regexp(txt,'</strong>','split');
            txt=cell2mat(txt);

            fprintf(fileID,'%s\n', txt);


            %Do the ranksum/t-test
            fprintf(1, ['\n\nRanksum or t-test p values for peak tPRP for '  handles.electrode_label{handles.peakLFPNo(ii_electrode)}  '\n']);
            fprintf(fileID, ['\n\nRanksum or t-test p values for peak tPRP for ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}  '\n']);

            [output_data] = drgMutiRanksumorTtest(input_data,fileID);

        end
    end


    %Now show tPRP peak for each electrode for laser window
    if these_plots_shown(6)==1
        for ii_electrode=1:4
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            hFig2 = figure(figNo);
            set(hFig2, 'units','normalized','position',[.1 .1 .3 .3])
            hold on

            bar_offset=0;

            glm_lfp_power=[];
            glm_ii_lfp=0;

            id_ii=0;
            input_data=[];

            for groupNo=1:length(groups)
                %Find the files for this mouse
                these_files=[];

                for fileNo=1:handles.no_files
                    if handles.group_no(fileNo)==groupNo
                        %Does this file exist?
                        this_filename=handles.FileName{fileNo};
                        FileName_out=['OBtoHIPtPRP_' this_filename(11:end)];
                        if exist([PathName_out FileName_out])~=0
                            these_files=[these_files fileNo];
                        end
                    end
                end

                %Get deltaP
                these_log_P_pre=[];
                these_log_P_laser=[];
                these_log_P_post=[];
                these_mice=[];

                if ~isempty(these_files)
                    for fileNo=these_files
                        %Plot the timecourse
                        % PathName_out=handles.PathName{fileNo};
                        this_filename=handles.FileName{fileNo};
                        FileName_out=['OBtoHIPtPRP_' this_filename(11:end)];
                        load([PathName_out FileName_out])
                        add_electrode=0;
                        if isempty(handles.file(fileNo).exclude_electrodes)
                            add_electrode=1;
                        else
                            if sum(handles.file(fileNo).exclude_electrodes==actual_electrode(ii_electrode))==0
                                add_electrode=1;
                            end

                        end
                        if add_electrode==1
                            log_P_timecourse=handles_out.electrode(ii_electrode).handles_out.mean_peakPower;
                            time=handles_out.electrode(ii_electrode).handles_out.time;
                            this_log_P_pre=mean(mean(log_P_timecourse(time<t_ref_end)));
                            these_log_P_pre=[these_log_P_pre this_log_P_pre];
                            this_log_P_laser=mean(mean(log_P_timecourse((time>=t_ref_end)&(time<=t_laser_end))));
                            these_log_P_laser=[these_log_P_laser this_log_P_laser];
                            this_log_P_post=mean(mean(log_P_timecourse(time>t_laser_end)));
                            these_log_P_post=[these_log_P_post this_log_P_post];
                            these_mice=[these_mice handles.mouse_no(fileNo)];
                        end
                        pffft=1;
                    end

                    these_log_P_pre=these_log_P_pre(~isnan(these_log_P_pre));
                    these_mice_pre=these_mice(~isnan(these_log_P_pre));
                    these_log_P_laser=these_log_P_laser(~isnan(these_log_P_laser));
                    these_mice_laser=these_mice(~isnan(these_log_P_laser));
                    these_log_P_post=these_log_P_post(~isnan(these_log_P_post));
                    these_mice_post=these_mice(~isnan(these_log_P_post));


                    %Laser
                    bar(bar_offset,mean(these_log_P_laser),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
                    [mean_out, CIout]=drgViolinPoint(these_log_P_laser,edges,bar_offset,rand_offset,'k','k',3);

                    glm_lfp_power.data(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=these_log_P_laser;
                    glm_lfp_power.time_window(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=1*ones(1,length(these_log_P_laser));
                    switch groupNo
                        case 1
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=1;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=1;
                        case 2
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=1;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=2;
                        case 3
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=2;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=1;
                        case 4
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=2;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=2;
                    end
                    glm_lfp_power.mouseNo(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=these_mice_laser;
                    glm_ii_lfp=glm_ii_lfp+length(these_log_P_laser);

                    id_ii=id_ii+1;
                    input_data(id_ii).data=these_log_P_laser;
                    input_data(id_ii).description=['laser ' group_label{groupNo}];

                  

                    bar_offset=bar_offset+1;
                else
                    bar_offset=bar_offset+3;
                end
            end

            title(['Peak tPRP for laser window ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}])


            xticks([0 1 2 3])
            xticklabels({group_label{1}, group_label{2}....
                ,group_label{3}, group_label{4}})


            % ylim([0 0.03])

            ylabel('dB')

            %Perform the glm for MI per mouse per odor pair
            fprintf(1, ['glm for peak tPRP for laser window ' handles.electrode_label{handles.peakLFPNo(ii_electrode)} '\n'])
            fprintf(fileID, ['glm for peak tPRP for laser window '  handles.electrode_label{handles.peakLFPNo(ii_electrode)} '\n']);

            % tbl = table(glm_lfp_power.data',glm_lfp_power.time_window',glm_lfp_power.group',...
            %     'VariableNames',{'logP','laser_post','group'});
            % mdl = fitglm(tbl,'logP~laser_post+group+laser_post*group'...
            %     ,'CategoricalVars',[2,3])

            % First, create a cell array of character labels corresponding to your group numbers
            % group_labels = {'Control', '5xFAD', 'Group3', 'Group4'}; % Add all your group labels here


            % for ii=1:length(re_arranged_groups)
            %     re_arranged_group_labels{ii}=handles.group_label{re_arranged_groups(ii)};
            % end
            tbl = table(glm_lfp_power.data',...
                categorical(glm_lfp_power.genotype', [1 2], genotype_label),...
                categorical(glm_lfp_power.laser', [1 2], laser_label),...
                'VariableNames',{'logP','genotype','laser'});
            mdl = fitglm(tbl,'logP~genotype+laser+genotype*laser')

            txt = evalc('mdl');
            txt=regexp(txt,'<strong>','split');
            txt=cell2mat(txt);
            txt=regexp(txt,'</strong>','split');
            txt=cell2mat(txt);

            fprintf(fileID,'%s\n', txt);


            %Do the ranksum/t-test
            fprintf(1, ['\n\nRanksum or t-test p values for peak tPRP for laser window for '  handles.electrode_label{handles.peakLFPNo(ii_electrode)}  '\n']);
            fprintf(fileID, ['\n\nRanksum or t-test p values for peak tPRP for laser window for ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}  '\n']);

            [output_data] = drgMutiRanksumorTtest(input_data,fileID);

        end
    end



    %Now show trough tPRP for each electrode
    if these_plots_shown(7)==1
        for ii_electrode=1:4
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            hFig2 = figure(figNo);
            set(hFig2, 'units','normalized','position',[.1 .1 .3 .3])
            hold on

            bar_offset=0;

            glm_lfp_power=[];
            glm_ii_lfp=0;

            id_ii=0;
            input_data=[];

            for groupNo=1:length(groups)
                %Find the files for this mouse
                these_files=[];

                for fileNo=1:handles.no_files
                    if handles.group_no(fileNo)==groupNo
                        %Does this file exist?
                        this_filename=handles.FileName{fileNo};
                        FileName_out=['OBtoHIPtPRP_' this_filename(11:end)];
                        if exist([PathName_out FileName_out])~=0
                            these_files=[these_files fileNo];
                        end
                    end
                end

                %Get deltaP
                these_log_P_pre=[];
                these_log_P_laser=[];
                these_log_P_post=[];
                these_mice=[];

                if ~isempty(these_files)
                    for fileNo=these_files
                        %Plot the timecourse
                        % PathName_out=handles.PathName{fileNo};
                        this_filename=handles.FileName{fileNo};
                        FileName_out=['OBtoHIPtPRP_' this_filename(11:end)];
                        load([PathName_out FileName_out])
                        add_electrode=0;
                        if isempty(handles.file(fileNo).exclude_electrodes)
                            add_electrode=1;
                        else
                            if sum(handles.file(fileNo).exclude_electrodes==actual_electrode(ii_electrode))==0
                                add_electrode=1;
                            end

                        end
                        if add_electrode==1
                            log_P_timecourse=handles_out.electrode(ii_electrode).handles_out.mean_troughPower;
                            time=handles_out.electrode(ii_electrode).handles_out.time;
                            this_log_P_pre=mean(mean(log_P_timecourse(time<t_ref_end)));
                            these_log_P_pre=[these_log_P_pre this_log_P_pre];
                            this_log_P_laser=mean(mean(log_P_timecourse((time>=t_ref_end)&(time<=t_laser_end))));
                            these_log_P_laser=[these_log_P_laser this_log_P_laser];
                            this_log_P_post=mean(mean(log_P_timecourse(time>t_laser_end)));
                            these_log_P_post=[these_log_P_post this_log_P_post];
                            these_mice=[these_mice handles.mouse_no(fileNo)];
                        end
                        pffft=1;
                    end

                    these_log_P_pre=these_log_P_pre(~isnan(these_log_P_pre));
                    these_mice_pre=these_mice(~isnan(these_log_P_pre));
                    these_log_P_laser=these_log_P_laser(~isnan(these_log_P_laser));
                    these_mice_laser=these_mice(~isnan(these_log_P_laser));
                    these_log_P_post=these_log_P_post(~isnan(these_log_P_post));
                    these_mice_post=these_mice(~isnan(these_log_P_post));


                    %Pre
                    bar(bar_offset,mean(these_log_P_pre),'LineWidth', 3,'EdgeColor','none','FaceColor',[31/255 31/255 140/255])
                    [mean_out, CIout]=drgViolinPoint(these_log_P_pre,edges,bar_offset,rand_offset,'k','k',3);

                    glm_lfp_power.data(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=these_log_P_pre;
                    glm_lfp_power.time_window(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=0*ones(1,length(these_log_P_pre));
                    switch groupNo
                        case 1
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=1;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=1;
                        case 2
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=1;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=2;
                        case 3
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=2;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=1;
                        case 4
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=2;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=2;
                    end
                    glm_lfp_power.mouseNo(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=these_mice_pre;
                    glm_ii_lfp=glm_ii_lfp+length(these_log_P_pre);

                    id_ii=id_ii+1;
                    input_data(id_ii).data=these_log_P_pre;
                    input_data(id_ii).description=['laser ' group_label{groupNo}];

                    %Laser
                    bar_offset=bar_offset+1;
                    bar(bar_offset,mean(these_log_P_laser),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
                    [mean_out, CIout]=drgViolinPoint(these_log_P_laser,edges,bar_offset,rand_offset,'k','k',3);

                    glm_lfp_power.data(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=these_log_P_laser;
                    glm_lfp_power.time_window(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=1*ones(1,length(these_log_P_laser));
                    switch groupNo
                        case 1
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=1;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=1;
                        case 2
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=1;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=2;
                        case 3
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=2;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=1;
                        case 4
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=2;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=2;
                    end
                    glm_lfp_power.mouseNo(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=these_mice_laser;
                    glm_ii_lfp=glm_ii_lfp+length(these_log_P_laser);

                    id_ii=id_ii+1;
                    input_data(id_ii).data=these_log_P_laser;
                    input_data(id_ii).description=['laser ' group_label{groupNo}];

                    %Post
                    bar_offset=bar_offset+1;
                    bar(bar_offset,mean(these_log_P_post),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
                    [mean_out, CIout]=drgViolinPoint(these_log_P_post,edges,bar_offset,rand_offset,'k','k',3);

                    glm_lfp_power.data(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=these_log_P_post;
                    glm_lfp_power.time_window(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=2*ones(1,length(these_log_P_post));
                    switch groupNo
                        case 1
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=1;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=1;
                        case 2
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=1;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=2;
                        case 3
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=2;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=1;
                        case 4
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=2;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=2;
                    end
                    glm_lfp_power.mouseNo(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=these_mice_post;
                    glm_ii_lfp=glm_ii_lfp+length(these_log_P_post);

                    id_ii=id_ii+1;
                    input_data(id_ii).data=these_log_P_post;
                    input_data(id_ii).description=['post ' group_label{groupNo}];

                    bar_offset=bar_offset+2;
                else
                    bar_offset=bar_offset+3;
                end
            end

            title(['Trough tPRP ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}])


            xticks([1 5 9 13])
            xticklabels({group_label{1}, group_label{2}....
                ,group_label{3}, group_label{4}})


            % ylim([0 0.03])

            ylabel('dB')

            %Perform the glm for MI per mouse per odor pair
            fprintf(1, ['glm for trough tPRP for ' handles.electrode_label{handles.peakLFPNo(ii_electrode)} '\n'])
            fprintf(fileID, ['glm for trough tPRP for '  handles.electrode_label{handles.peakLFPNo(ii_electrode)} '\n']);

            % tbl = table(glm_lfp_power.data',glm_lfp_power.time_window',glm_lfp_power.group',...
            %     'VariableNames',{'logP','laser_post','group'});
            % mdl = fitglm(tbl,'logP~laser_post+group+laser_post*group'...
            %     ,'CategoricalVars',[2,3])

            % First, create a cell array of character labels corresponding to your group numbers
            % group_labels = {'Control', '5xFAD', 'Group3', 'Group4'}; % Add all your group labels here


            % for ii=1:length(re_arranged_groups)
            %     re_arranged_group_labels{ii}=handles.group_label{re_arranged_groups(ii)};
            % end
            tbl = table(glm_lfp_power.data',...
                categorical(glm_lfp_power.time_window', [0 1 2], logP_time_window_labels),...
                categorical(glm_lfp_power.genotype', [1 2], genotype_label),...
                categorical(glm_lfp_power.laser', [1 2], laser_label),...
                'VariableNames',{'logP','time_window','genotype','laser'});
            mdl = fitglm(tbl,'logP~time_window+genotype+laser+time_window*genotype*laser')

            txt = evalc('mdl');
            txt=regexp(txt,'<strong>','split');
            txt=cell2mat(txt);
            txt=regexp(txt,'</strong>','split');
            txt=cell2mat(txt);

            fprintf(fileID,'%s\n', txt);


            %Do the ranksum/t-test
            fprintf(1, ['\n\nRanksum or t-test p values for trough tPRP for '  handles.electrode_label{handles.peakLFPNo(ii_electrode)}  '\n']);
            fprintf(fileID, ['\n\nRanksum or t-test p values for trough tPRP for ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}  '\n']);

            [output_data] = drgMutiRanksumorTtest(input_data,fileID);

        end
    end



    %Now show trough tPRP for each electrode for laser window only
    if these_plots_shown(8)==1
        for ii_electrode=1:4
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            hFig2 = figure(figNo);
            set(hFig2, 'units','normalized','position',[.1 .1 .3 .3])
            hold on

            bar_offset=0;

            glm_lfp_power=[];
            glm_ii_lfp=0;

            id_ii=0;
            input_data=[];

            for groupNo=1:length(groups)
                %Find the files for this mouse
                these_files=[];

                for fileNo=1:handles.no_files
                    if handles.group_no(fileNo)==groupNo
                        %Does this file exist?
                        this_filename=handles.FileName{fileNo};
                        FileName_out=['OBtoHIPtPRP_' this_filename(11:end)];
                        if exist([PathName_out FileName_out])~=0
                            these_files=[these_files fileNo];
                        end
                    end
                end

                %Get deltaP
                these_log_P_pre=[];
                these_log_P_laser=[];
                these_log_P_post=[];
                these_mice=[];

                if ~isempty(these_files)
                    for fileNo=these_files
                        %Plot the timecourse
                        % PathName_out=handles.PathName{fileNo};
                        this_filename=handles.FileName{fileNo};
                        FileName_out=['OBtoHIPtPRP_' this_filename(11:end)];
                        load([PathName_out FileName_out])
                        add_electrode=0;
                        if isempty(handles.file(fileNo).exclude_electrodes)
                            add_electrode=1;
                        else
                            if sum(handles.file(fileNo).exclude_electrodes==actual_electrode(ii_electrode))==0
                                add_electrode=1;
                            end

                        end
                        if add_electrode==1
                            log_P_timecourse=handles_out.electrode(ii_electrode).handles_out.mean_troughPower;
                            time=handles_out.electrode(ii_electrode).handles_out.time;
                            this_log_P_pre=mean(mean(log_P_timecourse(time<t_ref_end)));
                            these_log_P_pre=[these_log_P_pre this_log_P_pre];
                            this_log_P_laser=mean(mean(log_P_timecourse((time>=t_ref_end)&(time<=t_laser_end))));
                            these_log_P_laser=[these_log_P_laser this_log_P_laser];
                            this_log_P_post=mean(mean(log_P_timecourse(time>t_laser_end)));
                            these_log_P_post=[these_log_P_post this_log_P_post];
                            these_mice=[these_mice handles.mouse_no(fileNo)];
                        end
                        pffft=1;
                    end

                    these_log_P_pre=these_log_P_pre(~isnan(these_log_P_pre));
                    these_mice_pre=these_mice(~isnan(these_log_P_pre));
                    these_log_P_laser=these_log_P_laser(~isnan(these_log_P_laser));
                    these_mice_laser=these_mice(~isnan(these_log_P_laser));
                    these_log_P_post=these_log_P_post(~isnan(these_log_P_post));
                    these_mice_post=these_mice(~isnan(these_log_P_post));


                  

                    %Laser
       
                    bar(bar_offset,mean(these_log_P_laser),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
                    [mean_out, CIout]=drgViolinPoint(these_log_P_laser,edges,bar_offset,rand_offset,'k','k',3);

                    glm_lfp_power.data(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=these_log_P_laser;
                    glm_lfp_power.time_window(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=1*ones(1,length(these_log_P_laser));
                    switch groupNo
                        case 1
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=1;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=1;
                        case 2
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=1;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=2;
                        case 3
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=2;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=1;
                        case 4
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=2;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=2;
                    end
                    glm_lfp_power.mouseNo(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=these_mice_laser;
                    glm_ii_lfp=glm_ii_lfp+length(these_log_P_laser);

                    id_ii=id_ii+1;
                    input_data(id_ii).data=these_log_P_laser;
                    input_data(id_ii).description=['laser ' group_label{groupNo}];


                    bar_offset=bar_offset+1;
                   
                else
                    bar_offset=bar_offset+3;
                end
            end

            title(['Trough tPRP ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}])


            xticks([0 1 2 3])
            xticklabels({group_label{1}, group_label{2}....
                ,group_label{3}, group_label{4}})


            % ylim([0 0.03])

            ylabel('dB')

            %Perform the glm for MI per mouse per odor pair
            fprintf(1, ['glm for trough tPRP for laser window ' handles.electrode_label{handles.peakLFPNo(ii_electrode)} '\n'])
            fprintf(fileID, ['glm for trough tPRP for laser window '  handles.electrode_label{handles.peakLFPNo(ii_electrode)} '\n']);

            % tbl = table(glm_lfp_power.data',glm_lfp_power.time_window',glm_lfp_power.group',...
            %     'VariableNames',{'logP','laser_post','group'});
            % mdl = fitglm(tbl,'logP~laser_post+group+laser_post*group'...
            %     ,'CategoricalVars',[2,3])

            % First, create a cell array of character labels corresponding to your group numbers
            % group_labels = {'Control', '5xFAD', 'Group3', 'Group4'}; % Add all your group labels here


            % for ii=1:length(re_arranged_groups)
            %     re_arranged_group_labels{ii}=handles.group_label{re_arranged_groups(ii)};
            % end
            tbl = table(glm_lfp_power.data',...
                categorical(glm_lfp_power.genotype', [1 2], genotype_label),...
                categorical(glm_lfp_power.laser', [1 2], laser_label),...
                'VariableNames',{'logP','genotype','laser'});
            mdl = fitglm(tbl,'logP~genotype+laser+genotype*laser')

            txt = evalc('mdl');
            txt=regexp(txt,'<strong>','split');
            txt=cell2mat(txt);
            txt=regexp(txt,'</strong>','split');
            txt=cell2mat(txt);

            fprintf(fileID,'%s\n', txt);


            %Do the ranksum/t-test
            fprintf(1, ['\n\nRanksum or t-test p values for trough tPRP for laser window '  handles.electrode_label{handles.peakLFPNo(ii_electrode)}  '\n']);
            fprintf(fileID, ['\n\nRanksum or t-test p values for trough tPRP for laser window ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}  '\n']);

            [output_data] = drgMutiRanksumorTtest(input_data,fileID);

        end
    end

    % 
    % %Now show modulation index for each electrode pand each time window
    % edges=[0:0.01:1];
    % rand_offset=0.8;
    % for ii_electrode=1:4
    %     figNo=figNo+1;
    %     try
    %         close(figNo)
    %     catch
    %     end
    %     hFig2 = figure(figNo);
    %     set(hFig2, 'units','normalized','position',[.1 .1 .3 .3])
    %     hold on
    % 
    %     bar_offset=0;
    % 
    %     glm_lfp_power=[];
    %     glm_ii_lfp=0;
    % 
    %     id_ii=0;
    %     input_data=[];
    % 
    %     for groupNo=1:length(groups)
    %         %Find the files for this mouse
    %         these_files=[];
    %         these_days=[];
    % 
    %         for fileNo=1:handles.no_files
    %             if handles.group_no(fileNo)==groupNo
    %                 %Does this file exist?
    %                 this_filename=handles.FileName{fileNo};
    %                 FileName_out=['OBtoHIPtPRP_' this_filename(11:end)];
    %                 if exist([PathName_out FileName_out])~=0
    %                     these_files=[these_files fileNo];
    %                     these_days=[these_days handles.session_no(fileNo)];
    %                 end
    %             end
    %         end
    % 
    %         %Get deltaP
    %         these_log_P_pre=[];
    %         these_log_P_laser=[];
    %         these_log_P_post=[];
    %         these_mice=[];
    % 
    %         if ~isempty(these_files)
    %             for fileNo=these_files
    %                 %Plot the timecourse
    %                 % PathName_out=handles.PathName{fileNo};
    %                 this_filename=handles.FileName{fileNo};
    %                 FileName_out=['OBtoHIPtPRP_' this_filename(11:end)];
    %                 load([PathName_out FileName_out])
    %                 log_P_timecourse=handles_out.electrode(ii_electrode).handles_out.mod_indx;
    %                 time=handles_out.electrode(ii_electrode).handles_out.time;
    %                 this_log_P_pre=mean(mean(log_P_timecourse(time<t_ref_end)));
    %                 these_log_P_pre=[these_log_P_pre this_log_P_pre];
    %                 this_log_P_laser=mean(mean(log_P_timecourse((time>=t_ref_end)&(time<=t_laser_end))));
    %                 these_log_P_laser=[these_log_P_laser this_log_P_laser];
    %                 this_log_P_post=mean(mean(log_P_timecourse(time>t_laser_end)));
    %                 these_log_P_post=[these_log_P_post this_log_P_post];
    %                 these_mice=[these_mice handles.mouse_no(fileNo)];
    %                 pffft=1;
    %             end
    % 
    %             these_log_P_pre=these_log_P_pre(~isnan(these_log_P_pre));
    %             these_mice_pre=these_mice(~isnan(these_log_P_pre));
    %             these_log_P_laser=these_log_P_laser(~isnan(these_log_P_laser));
    %             these_mice_laser=these_mice(~isnan(these_log_P_laser));
    %             these_log_P_post=these_log_P_post(~isnan(these_log_P_post));
    %             these_mice_post=these_mice(~isnan(these_log_P_post));
    % 
    % 
    %             %Plot per day
    %             for ii_day=unique(these_days)
    %                 bar(bar_offset,mean(these_log_P_pre(these_days==ii_day)),'LineWidth', 3,'EdgeColor','none','FaceColor',[31/255 31/255 140/255])
    %                 [mean_out, CIout]=drgViolinPoint(these_log_P_pre,edges,bar_offset,rand_offset,'k','k',3);
    %                 bar_offset=bar_offset+1;
    % 
    %                 glm_lfp_power.data(glm_ii_lfp+1:glm_ii_lfp+sum(these_days==ii_day))=these_log_P_pre(these_days==ii_day);
    %                 glm_lfp_power.day(glm_ii_lfp+1:glm_ii_lfp+sum(these_days==ii_day))=ii_day*ones(1,sum(these_days==ii_day));
    %                 switch groupNo
    %                     case 1
    %                         glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+sum(these_days==ii_day))=1;
    %                         glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+sum(these_days==ii_day))=1;
    %                     case 2
    %                         glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+sum(these_days==ii_day))=1;
    %                         glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+sum(these_days==ii_day))=2;
    %                     case 3
    %                         glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+sum(these_days==ii_day))=2;
    %                         glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+sum(these_days==ii_day))=1;
    %                     case 4
    %                         glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+sum(these_days==ii_day))=2;
    %                         glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+sum(these_days==ii_day))=2;
    %                 end
    %                 glm_lfp_power.mouseNo(glm_ii_lfp+1:glm_ii_lfp+sum(these_days==ii_day))=these_mice_pre(these_days==ii_day);
    %                 glm_ii_lfp=glm_ii_lfp+sum(these_days==ii_day);
    % 
    %                 id_ii=id_ii+1;
    %                 input_data(id_ii).data=these_log_P_pre;
    %                 input_data(id_ii).description=['laser ' group_label{groupNo}];
    %             end
    % 
    % 
    % 
    %             bar_offset=bar_offset+2;
    %         else
    %             bar_offset=bar_offset+3;
    %         end
    %     end
    % 
    %     title(['Modulation index pre per day ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}])
    % 
    % 
    %     xticks([1 5 9 13])
    %     xticklabels({group_label{1}, group_label{2}....
    %         ,group_label{3}, group_label{4}})
    % 
    % 
    %     % ylim([0 0.03])
    % 
    %     ylabel('dB')
    % 
    %     %Perform the glm for MI per mouse per odor pair
    %     fprintf(1, ['glm for modulation index pre per day for ' handles.electrode_label{handles.peakLFPNo(ii_electrode)} '\n'])
    %     fprintf(fileID, ['glm for modulation index pre per day for '  handles.electrode_label{handles.peakLFPNo(ii_electrode)} '\n']);
    % 
    %     % tbl = table(glm_lfp_power.data',glm_lfp_power.time_window',glm_lfp_power.group',...
    %     %     'VariableNames',{'logP','laser_post','group'});
    %     % mdl = fitglm(tbl,'logP~laser_post+group+laser_post*group'...
    %     %     ,'CategoricalVars',[2,3])
    % 
    %     % First, create a cell array of character labels corresponding to your group numbers
    %     % group_labels = {'Control', '5xFAD', 'Group3', 'Group4'}; % Add all your group labels here
    % 
    % 
    %     % for ii=1:length(re_arranged_groups)
    %     %     re_arranged_group_labels{ii}=handles.group_label{re_arranged_groups(ii)};
    %     % end
    %     tbl = table(glm_lfp_power.data',glm_lfp_power.day',...
    %         categorical(glm_lfp_power.genotype', [1 2], genotype_label),...
    %         categorical(glm_lfp_power.laser', [1 2], laser_label),...
    %         'VariableNames',{'logP','day','genotype','laser'});
    %     mdl = fitglm(tbl,'logP~day+genotype+laser+day*genotype*laser')
    % 
    %     txt = evalc('mdl');
    %     txt=regexp(txt,'<strong>','split');
    %     txt=cell2mat(txt);
    %     txt=regexp(txt,'</strong>','split');
    %     txt=cell2mat(txt);
    % 
    %     fprintf(fileID,'%s\n', txt);
    % 
    % 
    %     %Do the ranksum/t-test
    %     fprintf(1, ['\n\nRanksum or t-test p values for modulation index for '  handles.electrode_label{handles.peakLFPNo(ii_electrode)}  '\n']);
    %     fprintf(fileID, ['\n\nRanksum or t-test p values for modulation index for ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}  '\n']);
    % 
    %     [output_data] = drgMutiRanksumorTtest(input_data,fileID);
    % 
    % end

    
    edges=[0:0.01:1];
    rand_offset=0.8;

    %Now show modulation index for each electrode for each time window
    if these_plots_shown(9)==1
        for ii_electrode=1:4
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            hFig2 = figure(figNo);
            set(hFig2, 'units','normalized','position',[.1 .1 .3 .3])
            hold on

            bar_offset=0;

            glm_lfp_power=[];
            glm_ii_lfp=0;

            id_ii=0;
            input_data=[];

            for groupNo=1:length(groups)
                %Find the files for this mouse
                these_files=[];

                for fileNo=1:handles.no_files
                    if handles.group_no(fileNo)==groupNo
                        %Does this file exist?
                        this_filename=handles.FileName{fileNo};
                        FileName_out=['OBtoHIPtPRP_' this_filename(11:end)];
                        if exist([PathName_out FileName_out])~=0
                            these_files=[these_files fileNo];
                        end
                    end
                end

                %Get deltaP
                these_log_P_pre=[];
                these_log_P_laser=[];
                these_log_P_post=[];
                these_mice=[];

                if ~isempty(these_files)
                    for fileNo=these_files
                        %Plot the timecourse
                        % PathName_out=handles.PathName{fileNo};
                        this_filename=handles.FileName{fileNo};
                        FileName_out=['OBtoHIPtPRP_' this_filename(11:end)];
                        load([PathName_out FileName_out])
                        add_electrode=0;
                        if isempty(handles.file(fileNo).exclude_electrodes)
                            add_electrode=1;
                        else
                            if sum(handles.file(fileNo).exclude_electrodes==actual_electrode(ii_electrode))==0
                                add_electrode=1;
                            end

                        end
                        if add_electrode==1
                            log_P_timecourse=handles_out.electrode(ii_electrode).handles_out.mod_indx;
                            time=handles_out.electrode(ii_electrode).handles_out.time;
                            this_log_P_pre=mean(mean(log_P_timecourse(time<t_ref_end)));
                            these_log_P_pre=[these_log_P_pre this_log_P_pre];
                            this_log_P_laser=mean(mean(log_P_timecourse((time>=t_ref_end)&(time<=t_laser_end))));
                            these_log_P_laser=[these_log_P_laser this_log_P_laser];
                            this_log_P_post=mean(mean(log_P_timecourse(time>t_laser_end)));
                            these_log_P_post=[these_log_P_post this_log_P_post];
                            these_mice=[these_mice handles.mouse_no(fileNo)];
                        end
                        pffft=1;
                    end

                    these_log_P_pre=these_log_P_pre(~isnan(these_log_P_pre));
                    these_mice_pre=these_mice(~isnan(these_log_P_pre));
                    these_log_P_laser=these_log_P_laser(~isnan(these_log_P_laser));
                    these_mice_laser=these_mice(~isnan(these_log_P_laser));
                    these_log_P_post=these_log_P_post(~isnan(these_log_P_post));
                    these_mice_post=these_mice(~isnan(these_log_P_post));


                    %Pre
                    bar(bar_offset,mean(these_log_P_pre),'LineWidth', 3,'EdgeColor','none','FaceColor',[31/255 31/255 140/255])
                    [mean_out, CIout]=drgViolinPoint(these_log_P_pre,edges,bar_offset,rand_offset,'k','k',3);

                    glm_lfp_power.data(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=these_log_P_pre;
                    glm_lfp_power.time_window(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=0*ones(1,length(these_log_P_pre));
                    switch groupNo
                        case 1
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=1;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=1;
                        case 2
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=1;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=2;
                        case 3
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=2;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=1;
                        case 4
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=2;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=2;
                    end
                    glm_lfp_power.mouseNo(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=these_mice_pre;
                    glm_ii_lfp=glm_ii_lfp+length(these_log_P_pre);

                    id_ii=id_ii+1;
                    input_data(id_ii).data=these_log_P_pre;
                    input_data(id_ii).description=['laser ' group_label{groupNo}];

                    %Laser
                    bar_offset=bar_offset+1;
                    bar(bar_offset,mean(these_log_P_laser),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
                    [mean_out, CIout]=drgViolinPoint(these_log_P_laser,edges,bar_offset,rand_offset,'k','k',3);

                    glm_lfp_power.data(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=these_log_P_laser;
                    glm_lfp_power.time_window(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=1*ones(1,length(these_log_P_laser));
                    switch groupNo
                        case 1
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=1;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=1;
                        case 2
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=1;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=2;
                        case 3
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=2;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=1;
                        case 4
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=2;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=2;
                    end
                    glm_lfp_power.mouseNo(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=these_mice_laser;
                    glm_ii_lfp=glm_ii_lfp+length(these_log_P_laser);

                    id_ii=id_ii+1;
                    input_data(id_ii).data=these_log_P_laser;
                    input_data(id_ii).description=['laser ' group_label{groupNo}];

                    %Post
                    bar_offset=bar_offset+1;
                    bar(bar_offset,mean(these_log_P_post),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
                    [mean_out, CIout]=drgViolinPoint(these_log_P_post,edges,bar_offset,rand_offset,'k','k',3);

                    glm_lfp_power.data(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=these_log_P_post;
                    glm_lfp_power.time_window(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=2*ones(1,length(these_log_P_post));
                    switch groupNo
                        case 1
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=1;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=1;
                        case 2
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=1;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=2;
                        case 3
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=2;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=1;
                        case 4
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=2;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=2;
                    end
                    glm_lfp_power.mouseNo(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=these_mice_post;
                    glm_ii_lfp=glm_ii_lfp+length(these_log_P_post);

                    id_ii=id_ii+1;
                    input_data(id_ii).data=these_log_P_post;
                    input_data(id_ii).description=['post ' group_label{groupNo}];

                    bar_offset=bar_offset+2;
                else
                    bar_offset=bar_offset+3;
                end
            end

            title(['Modulation index ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}])


            xticks([1 5 9 13])
            xticklabels({group_label{1}, group_label{2}....
                ,group_label{3}, group_label{4}})


            % ylim([0 0.03])

            ylabel('MI')

            %Perform the glm for MI per mouse per odor pair
            fprintf(1, ['glm for modulation index for ' handles.electrode_label{handles.peakLFPNo(ii_electrode)} '\n'])
            fprintf(fileID, ['glm for modulation index for '  handles.electrode_label{handles.peakLFPNo(ii_electrode)} '\n']);

            % tbl = table(glm_lfp_power.data',glm_lfp_power.time_window',glm_lfp_power.group',...
            %     'VariableNames',{'logP','laser_post','group'});
            % mdl = fitglm(tbl,'logP~laser_post+group+laser_post*group'...
            %     ,'CategoricalVars',[2,3])

            % First, create a cell array of character labels corresponding to your group numbers
            % group_labels = {'Control', '5xFAD', 'Group3', 'Group4'}; % Add all your group labels here


            % for ii=1:length(re_arranged_groups)
            %     re_arranged_group_labels{ii}=handles.group_label{re_arranged_groups(ii)};
            % end
            tbl = table(glm_lfp_power.data',...
                categorical(glm_lfp_power.time_window', [0 1 2], logP_time_window_labels),...
                categorical(glm_lfp_power.genotype', [1 2], genotype_label),...
                categorical(glm_lfp_power.laser', [1 2], laser_label),...
                'VariableNames',{'logP','time_window','genotype','laser'});
            mdl = fitglm(tbl,'logP~time_window+genotype+laser+time_window*genotype*laser')

            txt = evalc('mdl');
            txt=regexp(txt,'<strong>','split');
            txt=cell2mat(txt);
            txt=regexp(txt,'</strong>','split');
            txt=cell2mat(txt);

            fprintf(fileID,'%s\n', txt);


            %Do the ranksum/t-test
            fprintf(1, ['\n\nRanksum or t-test p values for modulation index for '  handles.electrode_label{handles.peakLFPNo(ii_electrode)}  '\n']);
            fprintf(fileID, ['\n\nRanksum or t-test p values for modulation index for ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}  '\n']);

            [output_data] = drgMutiRanksumorTtest(input_data,fileID);

        end
    end

    %Now show modulation index for each electrode for only one time window

    edges=[0:0.01:1];
    rand_offset=0.8;
    time_window=2; %1=Pre, 2=Laser, 3=Post
    if these_plots_shown(10)==1
        for ii_electrode=1:4
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            hFig2 = figure(figNo);
            set(hFig2, 'units','normalized','position',[.1 .1 .3 .3])
            hold on

            bar_offset=0;

            glm_lfp_power=[];
            glm_ii_lfp=0;

            id_ii=0;
            input_data=[];

            for groupNo=1:length(groups)
                %Find the files for this mouse
                these_files=[];

                for fileNo=1:handles.no_files
                    if handles.group_no(fileNo)==groupNo
                        %Does this file exist?
                        this_filename=handles.FileName{fileNo};
                        FileName_out=['OBtoHIPtPRP_' this_filename(11:end)];
                        if exist([PathName_out FileName_out])~=0
                            these_files=[these_files fileNo];
                        end
                    end
                end

                %Get deltaP
                these_log_P_pre=[];
                these_log_P_laser=[];
                these_log_P_post=[];
                these_mice=[];

                if ~isempty(these_files)
                    for fileNo=these_files
                        %Plot the timecourse
                        % PathName_out=handles.PathName{fileNo};
                        this_filename=handles.FileName{fileNo};
                        FileName_out=['OBtoHIPtPRP_' this_filename(11:end)];
                        load([PathName_out FileName_out])
                        add_electrode=0;
                        if isempty(handles.file(fileNo).exclude_electrodes)
                            add_electrode=1;
                        else
                            if sum(handles.file(fileNo).exclude_electrodes==actual_electrode(ii_electrode))==0
                                add_electrode=1;
                            end

                        end
                        if add_electrode==1
                            log_P_timecourse=handles_out.electrode(ii_electrode).handles_out.mod_indx;
                            time=handles_out.electrode(ii_electrode).handles_out.time;
                            this_log_P_pre=mean(mean(log_P_timecourse(time<t_ref_end)));
                            these_log_P_pre=[these_log_P_pre this_log_P_pre];
                            this_log_P_laser=mean(mean(log_P_timecourse((time>=t_ref_end)&(time<=t_laser_end))));
                            these_log_P_laser=[these_log_P_laser this_log_P_laser];
                            this_log_P_post=mean(mean(log_P_timecourse(time>t_laser_end)));
                            these_log_P_post=[these_log_P_post this_log_P_post];
                            these_mice=[these_mice handles.mouse_no(fileNo)];
                        end
                        pffft=1;
                    end

                    these_log_P_pre=these_log_P_pre(~isnan(these_log_P_pre));
                    these_mice_pre=these_mice(~isnan(these_log_P_pre));
                    these_log_P_laser=these_log_P_laser(~isnan(these_log_P_laser));
                    these_mice_laser=these_mice(~isnan(these_log_P_laser));
                    these_log_P_post=these_log_P_post(~isnan(these_log_P_post));
                    these_mice_post=these_mice(~isnan(these_log_P_post));


                    switch time_window
                        case 1
                            %Pre
                            bar(bar_offset,mean(these_log_P_pre),'LineWidth', 3,'EdgeColor','none','FaceColor',[31/255 31/255 140/255])
                            [mean_out, CIout]=drgViolinPoint(these_log_P_pre,edges,bar_offset,rand_offset,'k','k',3);

                            glm_lfp_power.data(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=these_log_P_pre;
                            glm_lfp_power.time_window(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=0*ones(1,length(these_log_P_pre));
                            switch groupNo
                                case 1
                                    glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=1;
                                    glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=1;
                                case 2
                                    glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=1;
                                    glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=2;
                                case 3
                                    glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=2;
                                    glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=1;
                                case 4
                                    glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=2;
                                    glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=2;
                            end
                            glm_lfp_power.mouseNo(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=these_mice_pre;
                            glm_ii_lfp=glm_ii_lfp+length(these_log_P_pre);

                            id_ii=id_ii+1;
                            input_data(id_ii).data=these_log_P_pre;
                            input_data(id_ii).description=['laser ' group_label{groupNo}];

                        case 2
                            %Laser
                           
                            bar(bar_offset,mean(these_log_P_laser),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
                            [mean_out, CIout]=drgViolinPoint(these_log_P_laser,edges,bar_offset,rand_offset,'k','k',3);

                            glm_lfp_power.data(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=these_log_P_laser;
                            glm_lfp_power.time_window(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=1*ones(1,length(these_log_P_laser));
                            switch groupNo
                                case 1
                                    glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=1;
                                    glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=1;
                                case 2
                                    glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=1;
                                    glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=2;
                                case 3
                                    glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=2;
                                    glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=1;
                                case 4
                                    glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=2;
                                    glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=2;
                            end
                            glm_lfp_power.mouseNo(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=these_mice_laser;
                            glm_ii_lfp=glm_ii_lfp+length(these_log_P_laser);

                            id_ii=id_ii+1;
                            input_data(id_ii).data=these_log_P_laser;
                            input_data(id_ii).description=['laser ' group_label{groupNo}];

                        case 3
                            %Post
                            bar_offset=bar_offset+1;
                            bar(bar_offset,mean(these_log_P_post),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
                            [mean_out, CIout]=drgViolinPoint(these_log_P_post,edges,bar_offset,rand_offset,'k','k',3);

                            glm_lfp_power.data(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=these_log_P_post;
                            glm_lfp_power.time_window(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=2*ones(1,length(these_log_P_post));
                            switch groupNo
                                case 1
                                    glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=1;
                                    glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=1;
                                case 2
                                    glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=1;
                                    glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=2;
                                case 3
                                    glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=2;
                                    glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=1;
                                case 4
                                    glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=2;
                                    glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=2;
                            end
                            glm_lfp_power.mouseNo(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=these_mice_post;
                            glm_ii_lfp=glm_ii_lfp+length(these_log_P_post);

                            id_ii=id_ii+1;
                            input_data(id_ii).data=these_log_P_post;
                            input_data(id_ii).description=['post ' group_label{groupNo}];
                    end
                    bar_offset=bar_offset+1;
                    
                else
                    bar_offset=bar_offset+3;
                end
            end

            title(['Modulation index ' handles.electrode_label{handles.peakLFPNo(ii_electrode)} ' ' logP_time_window_labels{time_window} ' only'])


            xticks([0 1 2 3])
            xticklabels({group_label{1}, group_label{2}....
                ,group_label{3}, group_label{4}})


            % ylim([0 0.03])

            ylabel('MI')

            %Perform the glm for MI per mouse per odor pair
            fprintf(1, ['glm for modulation index for ' handles.electrode_label{handles.peakLFPNo(ii_electrode)} ' ' logP_time_window_labels{time_window} ' only' '\n'])
            fprintf(fileID, ['glm for modulation index for '  handles.electrode_label{handles.peakLFPNo(ii_electrode)} ' ' logP_time_window_labels{time_window} ' only' '\n']);

            % tbl = table(glm_lfp_power.data',glm_lfp_power.time_window',glm_lfp_power.group',...
            %     'VariableNames',{'logP','laser_post','group'});
            % mdl = fitglm(tbl,'logP~laser_post+group+laser_post*group'...
            %     ,'CategoricalVars',[2,3])

            % First, create a cell array of character labels corresponding to your group numbers
            % group_labels = {'Control', '5xFAD', 'Group3', 'Group4'}; % Add all your group labels here


            % for ii=1:length(re_arranged_groups)
            %     re_arranged_group_labels{ii}=handles.group_label{re_arranged_groups(ii)};
            % end
            tbl = table(glm_lfp_power.data',...
                categorical(glm_lfp_power.genotype', [1 2], genotype_label),...
                categorical(glm_lfp_power.laser', [1 2], laser_label),...
                'VariableNames',{'logP','genotype','laser'});
            mdl = fitglm(tbl,'logP~genotype+laser+genotype*laser')

            txt = evalc('mdl');
            txt=regexp(txt,'<strong>','split');
            txt=cell2mat(txt);
            txt=regexp(txt,'</strong>','split');
            txt=cell2mat(txt);

            fprintf(fileID,'%s\n', txt);


            %Do the ranksum/t-test
            fprintf(1, ['\n\nRanksum or t-test p values for modulation index for '  handles.electrode_label{handles.peakLFPNo(ii_electrode)} ' ' logP_time_window_labels{time_window} ' only' '\n']);
            fprintf(fileID, ['\n\nRanksum or t-test p values for modulation index for ' handles.electrode_label{handles.peakLFPNo(ii_electrode)} ' ' logP_time_window_labels{time_window} ' only' '\n']);

            [output_data] = drgMutiRanksumorTtest(input_data,fileID);

        end
    end

    %Now show modulation index as a function of session number for each electrode for only one time window
    if these_plots_shown(11)==1
        edges=[0:0.01:1];
        rand_offset=0.8;
        time_window=2; %1=Pre, 2=Laser, 3=Post
        for ii_electrode=1:4
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            hFig2 = figure(figNo);
            set(hFig2, 'units','normalized','position',[.1 .1 .3 .3])
            hold on

            bar_offset=0;

            glm_lfp_power=[];
            glm_ii_lfp=0;

            id_ii=0;
            input_data=[];

            for groupNo=1:length(groups)
                %Find the files for this mouse
                these_files=[];

                for fileNo=1:handles.no_files
                    if handles.group_no(fileNo)==groupNo
                        %Does this file exist?
                        this_filename=handles.FileName{fileNo};
                        FileName_out=['OBtoHIPtPRP_' this_filename(11:end)];
                        if exist([PathName_out FileName_out])~=0
                            these_files=[these_files fileNo];
                        end
                    end
                end

                %Get deltaP
                these_log_P_pre=[];
                these_log_P_laser=[];
                these_log_P_post=[];
                these_mice=[];
                these_sessions=[];

                if ~isempty(these_files)
                    for fileNo=these_files
                        %Plot the timecourse
                        % PathName_out=handles.PathName{fileNo};
                        this_filename=handles.FileName{fileNo};
                        FileName_out=['OBtoHIPtPRP_' this_filename(11:end)];
                        load([PathName_out FileName_out])
                        add_electrode=0;
                        if isempty(handles.file(fileNo).exclude_electrodes)
                            add_electrode=1;
                        else
                            if sum(handles.file(fileNo).exclude_electrodes==actual_electrode(ii_electrode))==0
                                add_electrode=1;
                            end

                        end
                        if add_electrode==1
                            log_P_timecourse=handles_out.electrode(ii_electrode).handles_out.mod_indx;
                            time=handles_out.electrode(ii_electrode).handles_out.time;
                            this_log_P_pre=mean(mean(log_P_timecourse(time<t_ref_end)));
                            these_log_P_pre=[these_log_P_pre this_log_P_pre];
                            this_log_P_laser=mean(mean(log_P_timecourse((time>=t_ref_end)&(time<=t_laser_end))));
                            these_log_P_laser=[these_log_P_laser this_log_P_laser];
                            this_log_P_post=mean(mean(log_P_timecourse(time>t_laser_end)));
                            these_log_P_post=[these_log_P_post this_log_P_post];
                            these_mice=[these_mice handles.mouse_no(fileNo)];
                            these_sessions=[these_sessions handles.session_no(fileNo)];
                        end
                        pffft=1;
                    end

                    these_log_P_pre=these_log_P_pre(~isnan(these_log_P_pre));
                    these_mice_pre=these_mice(~isnan(these_log_P_pre));
                    these_log_P_laser=these_log_P_laser(~isnan(these_log_P_laser));
                    these_mice_laser=these_mice(~isnan(these_log_P_laser));
                    these_log_P_post=these_log_P_post(~isnan(these_log_P_post));
                    these_mice_post=these_mice(~isnan(these_log_P_post));


                    switch time_window
                        case 1
                            %Pre



                            glm_lfp_power.data(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=these_log_P_pre;
                            glm_lfp_power.time_window(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=0*ones(1,length(these_log_P_pre));
                            glm_lfp_power.session(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=these_sessions;
                            switch groupNo
                                case 1
                                    glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=1;
                                    glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=1;
                                    plot(these_sessions, these_log_P_pre,'ob')
                                case 2
                                    glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=1;
                                    glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=2;
                                    plot(these_sessions, these_log_P_pre,'xb')
                                case 3
                                    glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=2;
                                    glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=1;
                                    plot(these_sessions, these_log_P_pre,'*b')
                                case 4
                                    glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=2;
                                    glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=2;
                                    plot(these_sessions, these_log_P_pre,'sb')
                            end
                            glm_lfp_power.mouseNo(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_pre))=these_mice_pre;
                            glm_ii_lfp=glm_ii_lfp+length(these_log_P_pre);

                            id_ii=id_ii+1;
                            input_data(id_ii).data=these_log_P_pre;
                            input_data(id_ii).description=['laser ' group_label{groupNo}];

                        case 2
                            %Laser


                            glm_lfp_power.data(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=these_log_P_laser;
                            glm_lfp_power.time_window(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=1*ones(1,length(these_log_P_laser));
                            glm_lfp_power.session(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=these_sessions;
                            switch groupNo
                                case 1
                                    glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=1;
                                    glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=1;
                                    plot(these_sessions, these_log_P_laser,'ob')
                                case 2
                                    glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=1;
                                    glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=2;
                                    plot(these_sessions, these_log_P_laser,'xb')
                                case 3
                                    glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=2;
                                    glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=1;
                                    plot(these_sessions, these_log_P_laser,'*b')
                                case 4
                                    glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=2;
                                    glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=2;
                                    plot(these_sessions, these_log_P_laser,'sb')
                            end
                            glm_lfp_power.mouseNo(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_laser))=these_mice_laser;
                            glm_ii_lfp=glm_ii_lfp+length(these_log_P_laser);

                            id_ii=id_ii+1;
                            input_data(id_ii).data=these_log_P_laser;
                            input_data(id_ii).description=['laser ' group_label{groupNo}];

                        case 3
                            %Post


                            glm_lfp_power.data(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=these_log_P_post;
                            glm_lfp_power.time_window(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=2*ones(1,length(these_log_P_post));
                            glm_lfp_power.session(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=these_sessions;
                            switch groupNo
                                case 1
                                    glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=1;
                                    glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=1;
                                    plot(these_sessions, these_log_P_post,'ob')
                                case 2
                                    glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=1;
                                    glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=2;
                                    plot(these_sessions, these_log_P_post,'xb')
                                case 3
                                    glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=2;
                                    glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=1;
                                    plot(these_sessions, these_log_P_post,'*b')
                                case 4
                                    glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=2;
                                    glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=2;
                                    plot(these_sessions, these_log_P_post,'sb')
                            end
                            glm_lfp_power.mouseNo(glm_ii_lfp+1:glm_ii_lfp+length(these_log_P_post))=these_mice_post;
                            glm_ii_lfp=glm_ii_lfp+length(these_log_P_post);

                            id_ii=id_ii+1;
                            input_data(id_ii).data=these_log_P_post;
                            input_data(id_ii).description=['post ' group_label{groupNo}];
                    end

                    bar_offset=bar_offset+2;
                else
                    bar_offset=bar_offset+3;
                end
            end

            title(['Modulation index vs session ' handles.electrode_label{handles.peakLFPNo(ii_electrode)} ' ' logP_time_window_labels{time_window} ' only'])
            legend(group_label{1},group_label{2},group_label{3},group_label{4})
            xlabel('Session')
            ylabel('MI')


            %Perform the glm for MI per mouse per odor pair
            fprintf(1, ['glm for modulation index vs session for ' handles.electrode_label{handles.peakLFPNo(ii_electrode)} ' ' logP_time_window_labels{time_window} ' only' '\n'])
            fprintf(fileID, ['glm for modulation index vs session for '  handles.electrode_label{handles.peakLFPNo(ii_electrode)} ' ' logP_time_window_labels{time_window} ' only' '\n']);

            % tbl = table(glm_lfp_power.data',glm_lfp_power.time_window',glm_lfp_power.group',...
            %     'VariableNames',{'logP','laser_post','group'});
            % mdl = fitglm(tbl,'logP~laser_post+group+laser_post*group'...
            %     ,'CategoricalVars',[2,3])

            % First, create a cell array of character labels corresponding to your group numbers
            % group_labels = {'Control', '5xFAD', 'Group3', 'Group4'}; % Add all your group labels here


            % for ii=1:length(re_arranged_groups)
            %     re_arranged_group_labels{ii}=handles.group_label{re_arranged_groups(ii)};
            % end
            tbl = table(glm_lfp_power.data',...
                categorical(glm_lfp_power.genotype', [1 2], genotype_label),...
                categorical(glm_lfp_power.laser', [1 2], laser_label),...
                glm_lfp_power.session',...
                'VariableNames',{'logP','genotype','laser','session'});
            mdl = fitglm(tbl,'logP~genotype+laser+session+genotype*laser*session')

            txt = evalc('mdl');
            txt=regexp(txt,'<strong>','split');
            txt=cell2mat(txt);
            txt=regexp(txt,'</strong>','split');
            txt=cell2mat(txt);

            fprintf(fileID,'%s\n', txt);


            %Do the ranksum/t-test
            fprintf(1, ['\n\nRanksum or t-test p values for modulation index vs session for '  handles.electrode_label{handles.peakLFPNo(ii_electrode)} ' ' logP_time_window_labels{time_window} ' only' '\n']);
            fprintf(fileID, ['\n\nRanksum or t-test p values for modulation index vs session for ' handles.electrode_label{handles.peakLFPNo(ii_electrode)} ' ' logP_time_window_labels{time_window} ' only' '\n']);

            [output_data] = drgMutiRanksumorTtest(input_data,fileID);

        end
    end

    

    edges=[0:0.01:1];
    rand_offset=0.8;

    %Now show bar graph for peak angle variance for all windows
    if these_plots_shown(12)==1
        for ii_electrode=1:4
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            hFig2 = figure(figNo);
            set(hFig2, 'units','normalized','position',[.1 .1 .3 .3])
            hold on

            bar_offset=0;

            glm_lfp_power=[];
            glm_ii_lfp=0;

            id_ii=0;
            input_data=[];

            for groupNo=1:length(groups)
                %Find the files for this mouse
                these_files=[];

                for fileNo=1:handles.no_files
                    if handles.group_no(fileNo)==groupNo
                        %Does this file exist?
                        this_filename=handles.FileName{fileNo};
                        FileName_out=['OBtoHIPtPRP_' this_filename(11:end)];
                        if exist([PathName_out FileName_out])~=0
                            these_files=[these_files fileNo];
                        end
                    end
                end

                %Get deltaP
                these_peakAngle_var_pre=[];
                these_peakAngle_var_laser=[];
                these_peakAngle_var_post=[];
                these_mice=[];

                if ~isempty(these_files)
                    for fileNo=these_files
                        %Plot the timecourse
                        % PathName_out=handles.PathName{fileNo};
                        this_filename=handles.FileName{fileNo};
                        FileName_out=['OBtoHIPtPRP_' this_filename(11:end)];
                        load([PathName_out FileName_out])
                        add_electrode=0;
                        if isempty(handles.file(fileNo).exclude_electrodes)
                            add_electrode=1;
                        else
                            if sum(handles.file(fileNo).exclude_electrodes==actual_electrode(ii_electrode))==0
                                add_electrode=1;
                            end

                        end
                        if add_electrode==1
                            these_peakAngles=handles_out.electrode(ii_electrode).handles_out.peakAngle;
                            time=handles_out.electrode(ii_electrode).handles_out.time;
                            this_peakAngle_var_pre=((180/pi)^2)*circ_var(these_peakAngles(time<t_ref_end)'*pi/180);
                            these_peakAngle_var_pre=[these_peakAngle_var_pre this_peakAngle_var_pre];
                            this_peakAngle_var_laser=((180/pi)^2)*circ_var(these_peakAngles((time>=t_ref_end)&(time<=t_laser_end))'*pi/180);
                            these_peakAngle_var_laser=[these_peakAngle_var_laser this_peakAngle_var_laser];
                            this_peakAngle_var_post=((180/pi)^2)*circ_var(these_peakAngles(time>t_laser_end)'*pi/180);
                            these_peakAngle_var_post=[these_peakAngle_var_post this_peakAngle_var_post];
                            these_mice=[these_mice handles.mouse_no(fileNo)];
                        end
                        pffft=1;
                    end

                    these_peakAngle_var_pre=these_peakAngle_var_pre(~isnan(these_peakAngle_var_pre));
                    these_mice_pre=these_mice(~isnan(these_peakAngle_var_pre));
                    these_peakAngle_var_laser=these_peakAngle_var_laser(~isnan(these_peakAngle_var_laser));
                    these_mice_laser=these_mice(~isnan(these_peakAngle_var_laser));
                    these_peakAngle_var_post=these_peakAngle_var_post(~isnan(these_peakAngle_var_post));
                    these_mice_post=these_mice(~isnan(these_peakAngle_var_post));


                    %Pre
                    bar(bar_offset,mean(these_peakAngle_var_pre),'LineWidth', 3,'EdgeColor','none','FaceColor',[31/255 31/255 140/255])
                    [mean_out, CIout]=drgViolinPoint(these_peakAngle_var_pre,edges,bar_offset,rand_offset,'k','k',3);

                    glm_lfp_power.data(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_pre))=these_peakAngle_var_pre;
                    glm_lfp_power.time_window(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_pre))=0*ones(1,length(these_peakAngle_var_pre));
                    switch groupNo
                        case 1
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_pre))=1;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_pre))=1;
                        case 2
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_pre))=1;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_pre))=2;
                        case 3
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_pre))=2;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_pre))=1;
                        case 4
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_pre))=2;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_pre))=2;
                    end
                    glm_lfp_power.mouseNo(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_pre))=these_mice_pre;
                    glm_ii_lfp=glm_ii_lfp+length(these_peakAngle_var_pre);

                    id_ii=id_ii+1;
                    input_data(id_ii).data=these_peakAngle_var_pre;
                    input_data(id_ii).description=['laser ' group_label{groupNo}];

                    %Laser
                    bar_offset=bar_offset+1;
                    bar(bar_offset,mean(these_peakAngle_var_laser),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
                    [mean_out, CIout]=drgViolinPoint(these_peakAngle_var_laser,edges,bar_offset,rand_offset,'k','k',3);

                    glm_lfp_power.data(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_laser))=these_peakAngle_var_laser;
                    glm_lfp_power.time_window(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_laser))=1*ones(1,length(these_peakAngle_var_laser));
                    switch groupNo
                        case 1
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_laser))=1;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_laser))=1;
                        case 2
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_laser))=1;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_laser))=2;
                        case 3
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_laser))=2;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_laser))=1;
                        case 4
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_laser))=2;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_laser))=2;
                    end
                    glm_lfp_power.mouseNo(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_laser))=these_mice_laser;
                    glm_ii_lfp=glm_ii_lfp+length(these_peakAngle_var_laser);

                    id_ii=id_ii+1;
                    input_data(id_ii).data=these_peakAngle_var_laser;
                    input_data(id_ii).description=['laser ' group_label{groupNo}];

                    %Post
                    bar_offset=bar_offset+1;
                    bar(bar_offset,mean(these_peakAngle_var_post),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
                    [mean_out, CIout]=drgViolinPoint(these_peakAngle_var_post,edges,bar_offset,rand_offset,'k','k',3);

                    glm_lfp_power.data(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_post))=these_peakAngle_var_post;
                    glm_lfp_power.time_window(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_post))=2*ones(1,length(these_peakAngle_var_post));
                    switch groupNo
                        case 1
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_post))=1;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_post))=1;
                        case 2
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_post))=1;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_post))=2;
                        case 3
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_post))=2;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_post))=1;
                        case 4
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_post))=2;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_post))=2;
                    end
                    glm_lfp_power.mouseNo(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_post))=these_mice_post;
                    glm_ii_lfp=glm_ii_lfp+length(these_peakAngle_var_post);

                    id_ii=id_ii+1;
                    input_data(id_ii).data=these_peakAngle_var_post;
                    input_data(id_ii).description=['post ' group_label{groupNo}];

                    bar_offset=bar_offset+2;
                else
                    bar_offset=bar_offset+3;
                end
            end

            title(['PA variance ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}])


            xticks([1 5 9 13])
            xticklabels({group_label{1}, group_label{2}....
                ,group_label{3}, group_label{4}})


            % ylim([0 0.03])

            ylabel('degrees^2')

            %Perform the glm for MI per mouse per odor pair
            fprintf(1, ['glm for PA variance for ' handles.electrode_label{handles.peakLFPNo(ii_electrode)} '\n'])
            fprintf(fileID, ['glm for PA variance for '  handles.electrode_label{handles.peakLFPNo(ii_electrode)} '\n']);

            % tbl = table(glm_lfp_power.data',glm_lfp_power.time_window',glm_lfp_power.group',...
            %     'VariableNames',{'logP','laser_post','group'});
            % mdl = fitglm(tbl,'logP~laser_post+group+laser_post*group'...
            %     ,'CategoricalVars',[2,3])

            % First, create a cell array of character labels corresponding to your group numbers
            % group_labels = {'Control', '5xFAD', 'Group3', 'Group4'}; % Add all your group labels here


            % for ii=1:length(re_arranged_groups)
            %     re_arranged_group_labels{ii}=handles.group_label{re_arranged_groups(ii)};
            % end
            tbl = table(glm_lfp_power.data',...
                categorical(glm_lfp_power.time_window', [0 1 2], logP_time_window_labels),...
                categorical(glm_lfp_power.genotype', [1 2], genotype_label),...
                categorical(glm_lfp_power.laser', [1 2], laser_label),...
                'VariableNames',{'logP','time_window','genotype','laser'});
            mdl = fitglm(tbl,'logP~time_window+genotype+laser+time_window*genotype*laser')

            txt = evalc('mdl');
            txt=regexp(txt,'<strong>','split');
            txt=cell2mat(txt);
            txt=regexp(txt,'</strong>','split');
            txt=cell2mat(txt);

            fprintf(fileID,'%s\n', txt);


            %Do the ranksum/t-test
            fprintf(1, ['\n\nRanksum or t-test p values for PA variance for '  handles.electrode_label{handles.peakLFPNo(ii_electrode)}  '\n']);
            fprintf(fileID, ['\n\nRanksum or t-test p values for PA variance for ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}  '\n']);

            [output_data] = drgMutiRanksumorTtest(input_data,fileID);

        end
    end

    %Now show bar graph for peak angle variance for laser window
    if these_plots_shown(13)==1
        for ii_electrode=1:4
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            hFig2 = figure(figNo);
            set(hFig2, 'units','normalized','position',[.1 .1 .3 .3])
            hold on

            bar_offset=0;

            glm_lfp_power=[];
            glm_ii_lfp=0;

            id_ii=0;
            input_data=[];

            for groupNo=1:length(groups)
                %Find the files for this mouse
                these_files=[];

                for fileNo=1:handles.no_files
                    if handles.group_no(fileNo)==groupNo
                        %Does this file exist?
                        this_filename=handles.FileName{fileNo};
                        FileName_out=['OBtoHIPtPRP_' this_filename(11:end)];
                        if exist([PathName_out FileName_out])~=0
                            these_files=[these_files fileNo];
                        end
                    end
                end

                %Get deltaP
                these_peakAngle_var_pre=[];
                these_peakAngle_var_laser=[];
                these_peakAngle_var_post=[];
                these_mice=[];

                if ~isempty(these_files)
                    for fileNo=these_files
                        %Plot the timecourse
                        % PathName_out=handles.PathName{fileNo};
                        this_filename=handles.FileName{fileNo};
                        FileName_out=['OBtoHIPtPRP_' this_filename(11:end)];
                        load([PathName_out FileName_out])
                        add_electrode=0;
                        if isempty(handles.file(fileNo).exclude_electrodes)
                            add_electrode=1;
                        else
                            if sum(handles.file(fileNo).exclude_electrodes==actual_electrode(ii_electrode))==0
                                add_electrode=1;
                            end

                        end
                        if add_electrode==1
                            these_peakAngles=handles_out.electrode(ii_electrode).handles_out.peakAngle;
                            time=handles_out.electrode(ii_electrode).handles_out.time;
                            this_peakAngle_var_pre=((180/pi)^2)*circ_var(these_peakAngles(time<t_ref_end)'*pi/180);
                            these_peakAngle_var_pre=[these_peakAngle_var_pre this_peakAngle_var_pre];
                            this_peakAngle_var_laser=((180/pi)^2)*circ_var(these_peakAngles((time>=t_ref_end)&(time<=t_laser_end))'*pi/180);
                            these_peakAngle_var_laser=[these_peakAngle_var_laser this_peakAngle_var_laser];
                            this_peakAngle_var_post=((180/pi)^2)*circ_var(these_peakAngles(time>t_laser_end)'*pi/180);
                            these_peakAngle_var_post=[these_peakAngle_var_post this_peakAngle_var_post];
                            these_mice=[these_mice handles.mouse_no(fileNo)];
                        end
                        pffft=1;
                    end

                    these_peakAngle_var_pre=these_peakAngle_var_pre(~isnan(these_peakAngle_var_pre));
                    these_mice_pre=these_mice(~isnan(these_peakAngle_var_pre));
                    these_peakAngle_var_laser=these_peakAngle_var_laser(~isnan(these_peakAngle_var_laser));
                    these_mice_laser=these_mice(~isnan(these_peakAngle_var_laser));
                    these_peakAngle_var_post=these_peakAngle_var_post(~isnan(these_peakAngle_var_post));
                    these_mice_post=these_mice(~isnan(these_peakAngle_var_post));



                    %Laser
                    
                    bar(bar_offset,mean(these_peakAngle_var_laser),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
                    [mean_out, CIout]=drgViolinPoint(these_peakAngle_var_laser,edges,bar_offset,rand_offset,'k','k',3);

                    glm_lfp_power.data(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_laser))=these_peakAngle_var_laser;
                    glm_lfp_power.time_window(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_laser))=1*ones(1,length(these_peakAngle_var_laser));
                    switch groupNo
                        case 1
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_laser))=1;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_laser))=1;
                        case 2
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_laser))=1;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_laser))=2;
                        case 3
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_laser))=2;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_laser))=1;
                        case 4
                            glm_lfp_power.genotype(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_laser))=2;
                            glm_lfp_power.laser(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_laser))=2;
                    end
                    glm_lfp_power.mouseNo(glm_ii_lfp+1:glm_ii_lfp+length(these_peakAngle_var_laser))=these_mice_laser;
                    glm_ii_lfp=glm_ii_lfp+length(these_peakAngle_var_laser);

                    id_ii=id_ii+1;
                    input_data(id_ii).data=these_peakAngle_var_laser;
                    input_data(id_ii).description=['laser ' group_label{groupNo}];

                    bar_offset=bar_offset+1;

                else
                    bar_offset=bar_offset+3;
                end
            end

            title(['PA variance laser window ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}])


            xticks([0 1 2 3])
            xticklabels({group_label{1}, group_label{2}....
                ,group_label{3}, group_label{4}})


            % ylim([0 0.03])

            ylabel('degrees^2')

            %Perform the glm for MI per mouse per odor pair
            fprintf(1, ['glm for PA variance for laser window ' handles.electrode_label{handles.peakLFPNo(ii_electrode)} '\n'])
            fprintf(fileID, ['glm for PA variance for laser window '  handles.electrode_label{handles.peakLFPNo(ii_electrode)} '\n']);

            % tbl = table(glm_lfp_power.data',glm_lfp_power.time_window',glm_lfp_power.group',...
            %     'VariableNames',{'logP','laser_post','group'});
            % mdl = fitglm(tbl,'logP~laser_post+group+laser_post*group'...
            %     ,'CategoricalVars',[2,3])

            % First, create a cell array of character labels corresponding to your group numbers
            % group_labels = {'Control', '5xFAD', 'Group3', 'Group4'}; % Add all your group labels here


            % for ii=1:length(re_arranged_groups)
            %     re_arranged_group_labels{ii}=handles.group_label{re_arranged_groups(ii)};
            % end
            tbl = table(glm_lfp_power.data',...
                categorical(glm_lfp_power.genotype', [1 2], genotype_label),...
                categorical(glm_lfp_power.laser', [1 2], laser_label),...
                'VariableNames',{'logP','genotype','laser'});
            mdl = fitglm(tbl,'logP~genotype+laser+genotype*laser')

            txt = evalc('mdl');
            txt=regexp(txt,'<strong>','split');
            txt=cell2mat(txt);
            txt=regexp(txt,'</strong>','split');
            txt=cell2mat(txt);

            fprintf(fileID,'%s\n', txt);


            %Do the ranksum/t-test
            fprintf(1, ['\n\nRanksum or t-test p values for PA variance for laser window '  handles.electrode_label{handles.peakLFPNo(ii_electrode)}  '\n']);
            fprintf(fileID, ['\n\nRanksum or t-test p values for PA variance for laser window ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}  '\n']);

            [output_data] = drgMutiRanksumorTtest(input_data,fileID);

        end
    end


end


fclose(fileID);

pffft=1;





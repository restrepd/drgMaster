function drgOBtoHIP_summary_batch_LFP(choiceBatchPathName,choiceFileName)
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
t_ref_end=5; %In minutes
t_laser_end=65;


%This rearangement results in showing the groups in this order
%1 WT Untreated
%2 WT Treated
%3 5xFAD Untreated
%4 5xFAD Treated


if nargin==0
    [choiceFileName,choiceBatchPathName] = uigetfile({'drgbChoicesOBtoHIP*.m'},'Select the .m file with all the choices for analysis');
end

fprintf(1, ['\ndrgOBtoHIP_batch_LFP run for ' choiceFileName '\n\n']);

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
    FileName_out=['OBtoHIPLFP_' this_filename(11:end)]
    if exist([PathName_out FileName_out])==0
        fprintf(1, ['Program will be terminated because file No %d, does not exist\n'],filNum);
        all_files_present=0;
    end

end

%Text file for statistical output
fileID = fopen([choiceBatchPathName 'OBtoHIP_summary_batch_LFP_stats.txt'],'w');

tic


if all_files_present==1


    %Show delta logP bar graphs for each electrode
    these_mice=unique(handles.mouse_no);
    figNo=0;
    groups=unique(handles.group_no);
    edges=[-30:5:15];
    rand_offset=0.8;
    for ii_bw=1:length(LFPPowerFreqLow)

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
                        FileName_out=['OBtoHIPLFP_' this_filename(11:end)];
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
                        FileName_out=['OBtoHIPLFP_' this_filename(11:end)];
                        load([PathName_out FileName_out])
                        log_P_timecourse=handles_out.electrode(ii_electrode).handles_out.log_P_timecourse;
                        f=handles_out.electrode(ii_electrode).handles_out.f;
                        time=handles_out.electrode(ii_electrode).handles_out.time;
                        this_delta_log_P_laser=mean(mean(log_P_timecourse((f>=LFPPowerFreqLow(ii_bw))&(f<=LFPPowerFreqHigh(ii_bw)),(time>=t_ref_end)&(time<=t_laser_end))))-...
                            mean(mean(log_P_timecourse((f>=LFPPowerFreqLow(ii_bw))&(f<=LFPPowerFreqHigh(ii_bw)),time<=t_ref_end)));
                        these_delta_log_P_laser=[these_delta_log_P_laser this_delta_log_P_laser];
                        this_delta_log_P_post=mean(mean(log_P_timecourse((f>=LFPPowerFreqLow(ii_bw))&(f<=LFPPowerFreqHigh(ii_bw)),time>=t_laser_end)))-...
                            mean(mean(log_P_timecourse((f>=LFPPowerFreqLow(ii_bw))&(f<=LFPPowerFreqHigh(ii_bw)),time<=t_ref_end)));
                        these_delta_log_P_post=[these_delta_log_P_post this_delta_log_P_post];
                        these_mice=[these_mice handles.mouse_no(fileNo)];
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

            title(['LFP delta logP ' badwidth_names{ii_bw} ' ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}])


            xticks([0.5 3.5 6.5 9.5])
            xticklabels({group_label{1}, group_label{2}....
                ,group_label{3}, group_label{4}})


            % ylim([0 0.03])

            ylabel('dB')

            %Perform the glm for MI per mouse per odor pair
            fprintf(1, ['glm for delta logP for ' badwidth_names{ii_bw} ' ' handles.electrode_label{handles.peakLFPNo(ii_electrode)} '\n'])
            fprintf(fileID, ['glm for delta logP for ' badwidth_names{ii_bw} ' ' handles.electrode_label{handles.peakLFPNo(ii_electrode)} '\n']);

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
            fprintf(1, ['\n\nRanksum or t-test p values for delta logP for ' badwidth_names{ii_bw} ' ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}  '\n']);
            fprintf(fileID, ['\n\nRanksum or t-test p values for delta logP for ' badwidth_names{ii_bw} ' ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}  '\n']);

            [output_data] = drgMutiRanksumorTtest(input_data,fileID);

        end
    end

    %Now show log P for all bandwidths for each electrode
    these_mice=unique(handles.mouse_no);
   
    edges=[-30:5:15];
    rand_offset=0.8;
    for ii_bw=1:length(LFPPowerFreqLow)

        for ii_electrode=1:length(handles.peakLFPNo)
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

            for groupNo=1:length(unique(handles.group_no))
                %Find the files for this mouse
                these_files=[];

                 for fileNo=1:handles.no_files
                    if handles.group_no(fileNo)==groupNo
                        %Does this file exist?
                        this_filename=handles.FileName{fileNo};
                        FileName_out=['OBtoHIPLFP_' this_filename(11:end)];
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
                        FileName_out=['OBtoHIPLFP_' this_filename(11:end)];
                        load([PathName_out FileName_out])
                        log_P_timecourse=handles_out.electrode(ii_electrode).handles_out.log_P_timecourse;
                        f=handles_out.electrode(ii_electrode).handles_out.f;
                        time=handles_out.electrode(ii_electrode).handles_out.time;
                        this_log_P_pre=mean(mean(log_P_timecourse((f>=LFPPowerFreqLow(ii_bw))&(f<=LFPPowerFreqHigh(ii_bw)),time<t_ref_end)));
                        these_log_P_pre=[these_log_P_pre this_log_P_pre];
                        this_log_P_laser=mean(mean(log_P_timecourse((f>=LFPPowerFreqLow(ii_bw))&(f<=LFPPowerFreqHigh(ii_bw)),(time>=t_ref_end)&(time<=t_laser_end))));
                        these_log_P_laser=[these_log_P_laser this_log_P_laser];
                        this_log_P_post=mean(mean(log_P_timecourse((f>=LFPPowerFreqLow(ii_bw))&(f<=LFPPowerFreqHigh(ii_bw)),time>t_laser_end)));
                        these_log_P_post=[these_log_P_post this_log_P_post];
                        these_mice=[these_mice handles.mouse_no(fileNo)];
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

            title(['LFP logP ' badwidth_names{ii_bw} ' ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}])


            xticks([0.5 3.5 6.5 9.5])
            xticklabels({group_label{1}, group_label{2}....
                ,group_label{3}, group_label{4}})


            % ylim([0 0.03])

            ylabel('dB')

            %Perform the glm for MI per mouse per odor pair
            fprintf(1, ['glm for logP for ' badwidth_names{ii_bw} ' ' handles.electrode_label{handles.peakLFPNo(ii_electrode)} '\n'])
            fprintf(fileID, ['glm for logP for ' badwidth_names{ii_bw} ' ' handles.electrode_label{handles.peakLFPNo(ii_electrode)} '\n']);

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
            fprintf(1, ['\n\nRanksum or t-test p values for logP for ' badwidth_names{ii_bw} ' ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}  '\n']);
            fprintf(fileID, ['\n\nRanksum or t-test p values for logP for ' badwidth_names{ii_bw} ' ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}  '\n']);

            [output_data] = drgMutiRanksumorTtest(input_data,fileID);

        end
    end




    %Now show peak 40 Hz log P for each electrode
    these_mice=unique(handles.mouse_no);
    edges=[-30:5:15];
    rand_offset=0.8;

    show_f_bandwidth=[35 45];
    before_f=[22 27];
    after_f=[50 57];

    for ii_electrode=1:length(handles.peakLFPNo)
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
                    FileName_out=['OBtoHIPLFP_' this_filename(11:end)];
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
                    FileName_out=['OBtoHIPLFP_' this_filename(11:end)];
                    load([PathName_out FileName_out])
                    log_P_timecourse=handles_out.electrode(ii_electrode).handles_out.log_P_timecourse;
                    f=handles_out.electrode(ii_electrode).handles_out.f;
                    time=handles_out.electrode(ii_electrode).handles_out.time;

                    this_log_P_timecourse=zeros(1,size(log_P_timecourse,2));
                    this_log_P_timecourse(1,:)=mean(log_P_timecourse((f>=show_f_bandwidth(1))&(f<=show_f_bandwidth(2)),:));
                    f_b=mean(f((f>=show_f_bandwidth(1))&(f<=show_f_bandwidth(2))));
                    this_log_P_bef_timecourse=zeros(1,size(log_P_timecourse,2));
                    this_log_P_bef_timecourse(1,:)=mean(log_P_timecourse((f>=before_f(1))&(f<=before_f(2)),:));
                    bef_f=mean(f((f>=before_f(1))&(f<=before_f(2))));
                    this_log_P_af_timecourse=zeros(1,size(log_P_timecourse,2));
                    this_log_P_af_timecourse(1,:)=mean(log_P_timecourse((f>=after_f(1))&(f<=after_f(2)),:));
                    af_f=mean(f((f>=after_f(1))&(f<=after_f(2))));


                    slope=(this_log_P_af_timecourse-this_log_P_bef_timecourse)/(af_f-bef_f);
                    this_peak_log_P_af_timecourse=zeros(1,size(log_P_timecourse,2));

                    this_peak_log_P_af_timecourse=this_log_P_timecourse-(this_log_P_bef_timecourse+slope*(f_b-bef_f));

                    this_log_P_pre=mean(this_peak_log_P_af_timecourse(time<t_ref_end));
                    these_log_P_pre=[these_log_P_pre this_log_P_pre];
                    this_log_P_laser=mean(this_peak_log_P_af_timecourse((time>=t_ref_end)&(time<=t_laser_end)));
                    these_log_P_laser=[these_log_P_laser this_log_P_laser];
                    this_log_P_post=mean(this_peak_log_P_af_timecourse(time>t_laser_end));
                    these_log_P_post=[these_log_P_post this_log_P_post];
                    these_mice=[these_mice handles.mouse_no(fileNo)];

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

        title(['LFP peak 40 Hz logP ' badwidth_names{ii_bw} ' ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}])


        xticks([1 5 9 13])
        xticklabels({group_label{1}, group_label{2}....
                ,group_label{3}, group_label{4}})


        % ylim([0 0.03])

        ylabel('dB')

        %Perform the glm for MI per mouse per odor pair
        fprintf(1, ['glm for peak 40 Hz logP for '  handles.electrode_label{handles.peakLFPNo(ii_electrode)} '\n'])
        fprintf(fileID, ['glm for peak 40 Hz logP for '  handles.electrode_label{handles.peakLFPNo(ii_electrode)} '\n']);

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
        fprintf(1, ['\n\nRanksum or t-test p values for logP for ' badwidth_names{ii_bw} ' ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}  '\n']);
        fprintf(fileID, ['\n\nRanksum or t-test p values for logP for ' badwidth_names{ii_bw} ' ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}  '\n']);

        [output_data] = drgMutiRanksumorTtest(input_data,fileID);

    end

end



fclose(fileID);

pffft=1;





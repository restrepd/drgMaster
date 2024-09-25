function drgOBtoHIP_summary_batch_LFP(choiceBatchPathName,choiceFileName)
%Note: fitcnet will not work in Matlab versions earlier than 2021a

close all
clear all

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
    if iscell(handles.PathName)
        jtPathName=handles.PathName{filNum};
    else
        jtPathName=handles.PathName;
    end

    if exist([jtPathName jtFileName])==0
        fprintf(1, ['Program will be terminated because file No %d, does not exist\n'],filNum);
        all_files_present=0;
    end

end

%Text file for statistical output
fileID = fopen([choiceBatchPathName 'OBtoHIP_summary_batch_LFP_stats.txt'],'w');

tic

these_groups=[1 3 4 2];

if all_files_present==1


    %Show bar graphs for each electrode
    these_mice=unique(handles.mouse_no);
    figNo=0;
    groups=unique(handles.group_no);
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

            for groupNo=these_groups
                %Find the files for this mouse
                these_files=[];

                for fileNo=1:handles.no_files
                    if handles.group_no(fileNo)==groupNo
                        %Does this file exist?
                        PathName_out=handles.PathName{fileNo};
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
                        PathName_out=handles.PathName{fileNo};
                        this_filename=handles.FileName{fileNo};
                        FileName_out=['OBtoHIPLFP_' this_filename(11:end)];
                        load([PathName_out FileName_out])
                        log_P_timecourse=handles_out.electorde(ii_electrode).handles_out.log_P_timecourse;
                        f=handles_out.electorde(ii_electrode).handles_out.f;
                        this_delta_log_P_laser=mean(mean(log_P_timecourse((f>=LFPPowerFreqLow(ii_bw))&(f<=LFPPowerFreqHigh(ii_bw)),handles.ii_laser_start(fileNo):handles.ii_laser_end(fileNo))))-...
                            mean(mean(log_P_timecourse((f>=LFPPowerFreqLow(ii_bw))&(f<=LFPPowerFreqHigh(ii_bw)),1:handles.ii_laser_start(fileNo))));
                        these_delta_log_P_laser=[these_delta_log_P_laser this_delta_log_P_laser];
                        this_delta_log_P_post=mean(mean(log_P_timecourse((f>=LFPPowerFreqLow(ii_bw))&(f<=LFPPowerFreqHigh(ii_bw)),handles.ii_laser_end(fileNo)+1:end)))-...
                            mean(mean(log_P_timecourse((f>=LFPPowerFreqLow(ii_bw))&(f<=LFPPowerFreqHigh(ii_bw)),1:handles.ii_laser_start(fileNo))));
                        these_delta_log_P_post=[these_delta_log_P_post this_delta_log_P_post];
                        these_mice=[these_mice handles.mouse_no(fileNo)];
                        pffft=1;
                    end


                    %Laser
                    bar(bar_offset,mean(these_delta_log_P_laser),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
                    [mean_out, CIout]=drgViolinPoint(these_delta_log_P_laser,edges,bar_offset,rand_offset,'k','k',3);

                    glm_lfp_power.data(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=these_delta_log_P_laser;
                    glm_lfp_power.laser_post(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=0*ones(1,length(these_delta_log_P_laser));
                    glm_lfp_power.group(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=groupNo*ones(1,length(these_delta_log_P_laser));
                    glm_lfp_power.mouseNo(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_laser))=these_mice;
                    glm_ii_lfp=glm_ii_lfp+length(these_delta_log_P_laser);

                    id_ii=id_ii+1;
                    input_data(id_ii).data=these_delta_log_P_laser;
                    input_data(id_ii).description=['laser ' handles.group_label{groupNo}];

                    %Post
                    bar_offset=bar_offset+1;
                    bar(bar_offset,mean(these_delta_log_P_post),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
                    [mean_out, CIout]=drgViolinPoint(these_delta_log_P_post,edges,bar_offset,rand_offset,'k','k',3);

                    glm_lfp_power.data(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=these_delta_log_P_post;
                    glm_lfp_power.laser_post(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=0*ones(1,length(these_delta_log_P_post));
                    glm_lfp_power.group(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=groupNo*ones(1,length(these_delta_log_P_post));
                    glm_lfp_power.mouseNo(glm_ii_lfp+1:glm_ii_lfp+length(these_delta_log_P_post))=these_mice;
                    glm_ii_lfp=glm_ii_lfp+length(these_delta_log_P_post);

                    id_ii=id_ii+1;
                    input_data(id_ii).data=these_delta_log_P_post;
                    input_data(id_ii).description=['post ' handles.group_label{groupNo}];

                    bar_offset=bar_offset+2;
                else
                    bar_offset=bar_offset+3;
                end
            end

            title(['LFP delta logP ' badwidth_names{ii_bw} ' ' handles.electrode_label{handles.peakLFPNo(ii_electrode)}])


            xticks([0.5 3.5 6.5 9.5])
            xticklabels({handles.group_label{these_groups(1)}, handles.group_label{these_groups(2)},handles.group_label{these_groups(3)}, handles.group_label{these_groups(4)}})


            % ylim([0 0.03])

            ylabel('dB')

            %Perform the glm for MI per mouse per odor pair
            fprintf(1, ['glm for delta logP for ' badwidth_names{ii_bw} ' ' handles.electrode_label{handles.peakLFPNo(ii_electrode)} '\n'])
            fprintf(fileID, ['glm for delta logP for ' badwidth_names{ii_bw} ' ' handles.electrode_label{handles.peakLFPNo(ii_electrode)} '\n']);

            tbl = table(glm_lfp_power.data',glm_lfp_power.laser_post',glm_lfp_power.group',...
                'VariableNames',{'logP','laser_post','group'});
            mdl = fitglm(tbl,'logP~laser_post+group+laser_post*group'...
                ,'CategoricalVars',[2,3])

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

end


fclose(fileID);

pffft=1;





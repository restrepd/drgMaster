function drgAnalyzeBatchBehavior

%Ask user for the m file that contains information on what the user wants the analysis to be
%This file has all the information on what the user wants done, which files
%to process, what groups they fall into, etc
%
% An example of this file: drgbChoicesDanielPrelim
%
%

close all
clear all

[matFileName,matBatchPathName] = uigetfile({'*.mat'},'Select the .mat file with all the choices for analysis');
load([matBatchPathName matFileName])


%Plot percent correct
figNo = 1;

try
    close(figNo)
catch
end
hFig=figure(figNo);

set(hFig, 'units','normalized','position',[.02 .02 .95 .95])

max_session=max(handles.drgbchoices.session_no);
max_mouse=max(handles.drgbchoices.mouse_no);

for filNum=1:length(handles.drgbchoices.FileName)
    subplot(max_mouse,max_session,max_session*(handles.drgbchoices.mouse_no(filNum)-1)+handles.drgbchoices.session_no(filNum))
    set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
    % subplot(3,1,1)
    trials=1:length(handles.drgb.file(filNum).perCorr);
    
    %Plot in different colors
    plot(trials,handles.drgb.file(filNum).perCorr,'o','MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',2)
    hold on
    plot(trials(handles.drgb.file(filNum).encoding_trials),handles.drgb.file(filNum).perCorr(handles.drgb.file(filNum).encoding_trials),'o','MarkerEdgeColor',[0/255 158/255 115/255],'MarkerFaceColor',[0/255 158/255 115/255],'MarkerSize',2)
    plot(trials(handles.drgb.file(filNum).retrieval_trials),handles.drgb.file(filNum).perCorr(handles.drgb.file(filNum).retrieval_trials),'o','MarkerEdgeColor',[204/255 121/255 167/255],'MarkerFaceColor',[204/255 121/255 167/255],'MarkerSize',2)
     
    ylim([0 110]);
    try
        title([handles.drgbchoices.group_no_names{handles.drgbchoices.group_no(filNum)} ':' handles.drgbchoices.epoch_names{handles.drgbchoices.epoch(filNum)}])
    catch
        title([handles.drgbchoices.group_no_names{handles.drgbchoices.group_no(filNum)}])
    end
    remainder = rem((max_session*(handles.drgbchoices.mouse_no(filNum)-1)+handles.drgbchoices.session_no(filNum))-1,max_session);
    if  remainder==0
        ylabel(handles.drgbchoices.MouseName(handles.drgbchoices.mouse_no(filNum)))
    end
end

% title_str=inputdlg('Enter title');
%
% annotation('textbox', [0 0.9 1 0.1], ...
%     'String', title_str{1}, ...
%     'EdgeColor', 'none', ...
%     'HorizontalAlignment', 'center')

%Calculate the percent correct for the retreival trials

per_corr_per_mouse=zeros(max(handles.drgbchoices.mouse_no),4000);
no_trials_per_mouse=zeros(1,max(handles.drgbchoices.mouse_no));
group_no_per_mouse=zeros(1,max(handles.drgbchoices.mouse_no));
per_corr_distribution_per_mouse=zeros(max(handles.drgbchoices.mouse_no),length([0:5:100]));

group_legend{1}='WT';
group_legend{2}='Het';
group_legend{3}='KO';

per_corr_per_group=zeros(3,20000);
no_trials_per_group=zeros(1,3);



for filNum=1:length(handles.drgbchoices.FileName)
    
    trials=1:length(handles.drgb.file(filNum).perCorr);
    
    %     per_corr_per_mouse(handles.drgbchoices.mouse_no(filNum),no_trials_per_mouse(handles.drgbchoices.mouse_no(filNum))+1:no_trials_per_mouse(handles.drgbchoices.mouse_no(filNum))...
    %         +sum(handles.drgb.file(filNum).retrieval_trials))=handles.drgb.file(filNum).perCorr(handles.drgb.file(filNum).retrieval_trials);
    %
    %     no_trials_per_mouse(handles.drgbchoices.mouse_no(filNum))=no_trials_per_mouse(handles.drgbchoices.mouse_no(filNum))+sum(handles.drgb.file(filNum).retrieval_trials);
    %     group_no_per_mouse(handles.drgbchoices.mouse_no(filNum))=handles.drgbchoices.group_no_per_mouse(filNum);
    
    per_corr_per_mouse(handles.drgbchoices.mouse_no(filNum),no_trials_per_mouse(handles.drgbchoices.mouse_no(filNum))+1:no_trials_per_mouse(handles.drgbchoices.mouse_no(filNum))...
        +length(handles.drgb.file(filNum).perCorr))=handles.drgb.file(filNum).perCorr;
    
    no_trials_per_mouse(handles.drgbchoices.mouse_no(filNum))=no_trials_per_mouse(handles.drgbchoices.mouse_no(filNum))+length(handles.drgb.file(filNum).perCorr);
    group_no_per_mouse(handles.drgbchoices.mouse_no(filNum))=handles.drgbchoices.group_no(filNum);

    for pcc=0:length([0:5:100])-1 
        per_corr_distribution_per_mouse(handles.drgbchoices.mouse_no(filNum),pcc+1)=per_corr_distribution_per_mouse(handles.drgbchoices.mouse_no(filNum),pcc+1)+sum(handles.drgb.file(filNum).perCorr==5*pcc);
    end
    
    per_corr_per_group(handles.drgbchoices.group_no(filNum),no_trials_per_group(handles.drgbchoices.group_no(filNum))+1:no_trials_per_group(handles.drgbchoices.group_no(filNum))...
        +length(handles.drgb.file(filNum).perCorr))=handles.drgb.file(filNum).perCorr;
    
    no_trials_per_group(handles.drgbchoices.group_no(filNum))=no_trials_per_group(handles.drgbchoices.group_no(filNum))+length(handles.drgb.file(filNum).perCorr);
    
end

per_corr_per_mouse_per_mouse=zeros(1,max(handles.drgbchoices.mouse_no));

glm_pcorr=[];
glm_ii=0;

for ii=1:max(handles.drgbchoices.mouse_no)
    mean_per_corr_per_mouse(ii)=mean(per_corr_per_mouse(ii,1:no_trials_per_mouse(ii)));
    per_corr_distribution_per_mouse(ii,:)=per_corr_distribution_per_mouse(ii,:)/sum(per_corr_distribution_per_mouse(ii,:));
    glm_pcorr.data(glm_ii+1:glm_ii+length([0:5:100]))=per_corr_distribution_per_mouse(ii,:)/sum(per_corr_distribution_per_mouse(ii,:));
    glm_pcorr.per_ii(glm_ii+1:glm_ii+length([0:5:100]))=[0:5:100];
    glm_pcorr.grNo(glm_ii+1:glm_ii+length([0:5:100]))=group_no_per_mouse(ii);
    glm_ii=glm_ii+length([0:5:100]);
end

mean_per_corr_distribution_per_group=zeros(max(handles.drgbchoices.group_no),length([0:5:100]));
for grNo=1:max(handles.drgbchoices.group_no)
    mean_per_corr_distribution_per_group(grNo,:)=mean(per_corr_distribution_per_mouse(group_no_per_mouse==grNo,:),1);
end

figNo = figNo +1;

try
    close(figNo)
catch
end
hFig=figure(figNo);

set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
hold on


bar_offset = 0;
edges=[80:1:100];
rand_offset=0.8;
id_ii=0;
input_data=[];

for grNo=1:max(group_no_per_mouse)
    
    switch grNo
        case 1
            bar(bar_offset,mean(mean_per_corr_per_mouse(group_no_per_mouse==grNo)),'g','LineWidth', 3,'EdgeColor','none')
        case 2
            bar(bar_offset,mean(mean_per_corr_per_mouse(group_no_per_mouse==grNo)),'b','LineWidth', 3,'EdgeColor','none')
        case 3
            bar(bar_offset,mean(mean_per_corr_per_mouse(group_no_per_mouse==grNo)),'y','LineWidth', 3,'EdgeColor','none')
    end
    
     
    %Violin plot
    
    [mean_out, CIout]=drgViolinPoint(mean_per_corr_per_mouse(group_no_per_mouse==grNo)...
        ,edges,bar_offset,rand_offset,'k','k',1);
    
    %                                 %Save data for  ranksum
    
    id_ii=id_ii+1;
    input_data(id_ii).data=mean_per_corr_per_mouse(group_no_per_mouse==grNo);
    input_data(id_ii).description=group_legend{grNo};
    
    
    
    
    bar_offset = bar_offset + 1;
    
    
    
    
    
end

title(['Percent correct for retreival'])


xticks([1 2 3])
xticklabels({'WT', 'Het', 'KO'})


ylabel('Percent correct')



%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for percent correct\n'])
[output_data] = drgMutiRanksumorTtest(input_data);



figNo = figNo +1;

try
    close(figNo)
catch
end
hFig=figure(figNo);

set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
hold on


bar_offset = 0;
edges=[0:5:100];
rand_offset=0.8;
id_ii=0;
input_data=[];

for grNo=1:max(group_no_per_mouse)
    
    switch grNo
        case 1
            bar(bar_offset,mean(per_corr_per_group(grNo,1:no_trials_per_group(grNo))),'b','LineWidth', 3,'EdgeColor','none')
        case 2
            bar(bar_offset,mean(per_corr_per_group(grNo,1:no_trials_per_group(grNo))),'g','LineWidth', 3,'EdgeColor','none')
        case 3
            bar(bar_offset,mean(per_corr_per_group(grNo,1:no_trials_per_group(grNo))),'y','LineWidth', 3,'EdgeColor','none')
    end
    
     
    %Violin plot
    
    [mean_out, CIout]=drgViolinPoint(per_corr_per_group(grNo,1:no_trials_per_group(grNo))'...
        ,edges,bar_offset,rand_offset,'k','k',1);
    
    %                                 %Save data for  ranksum
    
    id_ii=id_ii+1;
    input_data(id_ii).data=per_corr_per_group(grNo,1:no_trials_per_group(grNo));
    input_data(id_ii).description=group_legend{grNo};
    
    
    
    
    bar_offset = bar_offset + 1;
    
    
    
    
    
end

title(['Percent correct per group'])


xticks([0 1 2])
xticklabels({'WT', 'Het', 'KO'})


ylabel('Percent correct')



%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for percent correct\n'])
[output_data] = drgMutiRanksumorTtest(input_data);

%Interesting, but is it real? Let's do per mouse distributionsfor grNo=1:max(group_no_per_mouse)

figNo = figNo +1;

try
    close(figNo)
catch
end
hFig=figure(figNo);

set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
hold on

for grNo=max(group_no_per_mouse):-1:1
    
    these_per_corr_distribution=[];
    these_per_corr_distribution=per_corr_distribution_per_mouse(group_no_per_mouse==grNo,:);
    
    this_mean_per_corr_distribution=mean_per_corr_distribution_per_group(grNo,:);
    
    CI=[];
    CI = bootci(1000, {@mean, these_per_corr_distribution})';
    CI(:,1)= this_mean_per_corr_distribution'-CI(:,1);
    CI(:,2)=CI(:,2)- this_mean_per_corr_distribution';
    
    
    switch grNo
        case 1
            [hlCR, hpCR] = boundedline([0:5:100],this_mean_per_corr_distribution', CI, 'b');
            %             plot([0:5:100], mean_per_corr_distribution_per_group(grNo,:), '-g','LineWidth', 3)
        case 2
            [hlCR, hpCR] = boundedline([0:5:100],this_mean_per_corr_distribution', CI, 'g');
        case 3
            [hlCR, hpCR] = boundedline([0:5:100],this_mean_per_corr_distribution', CI, 'y');
    end
end
                          
%Perform the glm
fprintf(1, ['\n\nglm for percent correct distribution \n'])
tbl = table(glm_pcorr.data',glm_pcorr.per_ii',glm_pcorr.grNo',...
    'VariableNames',{'fraction','pcorr','grNo'});
mdl = fitglm(tbl,'fraction~pcorr+grNo+pcorr*grNo','CategoricalVars',[3])

pfft=1;



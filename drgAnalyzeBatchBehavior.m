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

proficient_pCorr=80;
naive_pCorr=[45 65];
proficient_window=40;
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
itis_per_mouse=zeros(max(handles.drgbchoices.mouse_no),4000);
no_trials_per_mouse=zeros(1,max(handles.drgbchoices.mouse_no));
group_no_per_mouse=zeros(1,max(handles.drgbchoices.mouse_no));
per_corr_distribution_per_mouse=zeros(max(handles.drgbchoices.mouse_no),length([0:5:100]));

group_legend{1}='WT';
group_legend{2}='Het';
group_legend{3}='KO';

per_corr_per_group=zeros(3,20000);
no_trials_per_group=zeros(1,3);

no_sessions_to_proficient_per_mouse=zeros(1,max(handles.drgbchoices.mouse_no));
no_sessions_to_proficient_per_mouse_found=zeros(1,max(handles.drgbchoices.mouse_no));

for filNum=1:length(handles.drgbchoices.FileName)
    
    if length(handles.drgb.file(filNum).perCorr)>1
        %One of EAPA files has only one trial, exclude it
        trials=1:length(handles.drgb.file(filNum).perCorr);
        
        %     per_corr_per_mouse(handles.drgbchoices.mouse_no(filNum),no_trials_per_mouse(handles.drgbchoices.mouse_no(filNum))+1:no_trials_per_mouse(handles.drgbchoices.mouse_no(filNum))...
        %         +sum(handles.drgb.file(filNum).retrieval_trials))=handles.drgb.file(filNum).perCorr(handles.drgb.file(filNum).retrieval_trials);
        %
        %     no_trials_per_mouse(handles.drgbchoices.mouse_no(filNum))=no_trials_per_mouse(handles.drgbchoices.mouse_no(filNum))+sum(handles.drgb.file(filNum).retrieval_trials);
        %     group_no_per_mouse(handles.drgbchoices.mouse_no(filNum))=handles.drgbchoices.group_no_per_mouse(filNum);
        
        per_corr_per_mouse(handles.drgbchoices.mouse_no(filNum),no_trials_per_mouse(handles.drgbchoices.mouse_no(filNum))+1:no_trials_per_mouse(handles.drgbchoices.mouse_no(filNum))...
            +length(handles.drgb.file(filNum).perCorr))=handles.drgb.file(filNum).perCorr;
        
        if (length(handles.drgb.file(filNum).perCorr)>=proficient_window)&(no_sessions_to_proficient_per_mouse_found(1,handles.drgbchoices.mouse_no(filNum))==0)
            found_prof=0;
            ii=0;
            while (ii<length(handles.drgb.file(filNum).perCorr)-proficient_window)&(found_prof==0)
                ii=ii+1;
                if sum(handles.drgb.file(filNum).perCorr(ii:ii+proficient_window-1)>=proficient_pCorr)==proficient_window
                    found_prof=1;
                end
            end
            if found_prof==1
                no_sessions_to_proficient_per_mouse(1,handles.drgbchoices.mouse_no(filNum))=no_sessions_to_proficient_per_mouse(1,handles.drgbchoices.mouse_no(filNum))+(ii/length(handles.drgb.file(filNum).perCorr));
                no_sessions_to_proficient_per_mouse_found(1,handles.drgbchoices.mouse_no(filNum))=1;
            else
                no_sessions_to_proficient_per_mouse(1,handles.drgbchoices.mouse_no(filNum))=no_sessions_to_proficient_per_mouse(1,handles.drgbchoices.mouse_no(filNum))+1;
            end
        end
        
        
        these_itis=handles.drgb.file(filNum).times(2:end)-handles.drgb.file(filNum).times(1:end-1);
        %Note that the first iti is the same as the second
        these_itis=[these_itis(1) these_itis];
        itis_per_mouse(handles.drgbchoices.mouse_no(filNum),no_trials_per_mouse(handles.drgbchoices.mouse_no(filNum))+1:no_trials_per_mouse(handles.drgbchoices.mouse_no(filNum))...
            +length(handles.drgb.file(filNum).perCorr))=these_itis;
        
        
        no_trials_per_mouse(handles.drgbchoices.mouse_no(filNum))=no_trials_per_mouse(handles.drgbchoices.mouse_no(filNum))+length(handles.drgb.file(filNum).perCorr);
        group_no_per_mouse(handles.drgbchoices.mouse_no(filNum))=handles.drgbchoices.group_no(filNum);
        
        for pcc=0:length([0:5:100])-1
            per_corr_distribution_per_mouse(handles.drgbchoices.mouse_no(filNum),pcc+1)=per_corr_distribution_per_mouse(handles.drgbchoices.mouse_no(filNum),pcc+1)+sum(handles.drgb.file(filNum).perCorr==5*pcc);
        end
        
        per_corr_per_group(handles.drgbchoices.group_no(filNum),no_trials_per_group(handles.drgbchoices.group_no(filNum))+1:no_trials_per_group(handles.drgbchoices.group_no(filNum))...
            +length(handles.drgb.file(filNum).perCorr))=handles.drgb.file(filNum).perCorr;
        
        no_trials_per_group(handles.drgbchoices.group_no(filNum))=no_trials_per_group(handles.drgbchoices.group_no(filNum))+length(handles.drgb.file(filNum).perCorr);
    end
end

per_corr_per_mouse_per_mouse=zeros(1,max(handles.drgbchoices.mouse_no));

glm_pcorr=[];
glm_ii=0;

for ii=1:max(handles.drgbchoices.mouse_no)
    mean_per_corr_per_mouse(ii)=mean(per_corr_per_mouse(ii,1:no_trials_per_mouse(ii)));
    this_mouse_per_corr=per_corr_per_mouse(ii,1:no_trials_per_mouse(ii));
    mean_per_corr_per_mouse_prof(ii)=mean(this_mouse_per_corr(this_mouse_per_corr>=proficient_pCorr));
    fraction_per_corr_per_mouse_prof(ii)=sum(this_mouse_per_corr>=proficient_pCorr)/length(this_mouse_per_corr);
    ratio_per_corr_per_mouse_prof(ii)=sum(this_mouse_per_corr>=proficient_pCorr)/sum((this_mouse_per_corr<=naive_pCorr(2))&(this_mouse_per_corr>=naive_pCorr(1)));
    per_corr_distribution_per_mouse(ii,:)=per_corr_distribution_per_mouse(ii,:)/sum(per_corr_distribution_per_mouse(ii,:));
    glm_pcorr.data(glm_ii+1:glm_ii+length([0:5:100]))=per_corr_distribution_per_mouse(ii,:)/sum(per_corr_distribution_per_mouse(ii,:));
    glm_pcorr.per_ii(glm_ii+1:glm_ii+length([0:5:100]))=[0:5:100];
    glm_pcorr.grNo(glm_ii+1:glm_ii+length([0:5:100]))=group_no_per_mouse(ii);
    glm_ii=glm_ii+length([0:5:100]);
    mean_iti_per_mouse(ii)=mean(itis_per_mouse(ii,1:no_trials_per_mouse(ii)));
    this_mouse_itis=itis_per_mouse(ii,1:no_trials_per_mouse(ii));
    mean_iti_per_mouse_proficient(ii)=mean(this_mouse_itis(this_mouse_per_corr>=proficient_pCorr));
end

mean_per_corr_distribution_per_group=zeros(max(handles.drgbchoices.group_no),length([0:5:100]));
for grNo=1:max(handles.drgbchoices.group_no)
    mean_per_corr_distribution_per_group(grNo,:)=mean(per_corr_distribution_per_mouse(group_no_per_mouse==grNo,:),1);
end

%Generate mean percent correct per mouse bar graph
figNo = figNo +1;

try
    close(figNo)
catch
end
hFig=figure(figNo);

set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
hold on


bar_offset = 0;
edges=[proficient_pCorr:1:100];
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


xticks([0 1 2])
xticklabels({'WT', 'Het', 'KO'})


ylabel('Percent correct')



%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for percent correct\n'])
[output_data] = drgMutiRanksumorTtest(input_data);

%Generate mean proficient iti per mouse bar graph
id_ii=0;
input_data=[];

figNo = figNo +1;

try
    close(figNo)
catch
end
hFig=figure(figNo);

set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
hold on


bar_offset = 0;
edges=[0:2:50];
rand_offset=0.8;
id_ii=0;
input_data=[];

for grNo=1:max(group_no_per_mouse)
    
    switch grNo
        case 1
            bar(bar_offset,mean(mean_iti_per_mouse_proficient((group_no_per_mouse==grNo)&(~isnan(mean_iti_per_mouse_proficient)))),'g','LineWidth', 3,'EdgeColor','none')
        case 2
            bar(bar_offset,mean(mean_iti_per_mouse_proficient((group_no_per_mouse==grNo)&(~isnan(mean_iti_per_mouse_proficient)))),'b','LineWidth', 3,'EdgeColor','none')
        case 3
            bar(bar_offset,mean(mean_iti_per_mouse_proficient((group_no_per_mouse==grNo)&(~isnan(mean_iti_per_mouse_proficient)))),'y','LineWidth', 3,'EdgeColor','none')
    end
    
     
    %Violin plot
    
    [mean_out, CIout]=drgViolinPoint(mean_iti_per_mouse_proficient((group_no_per_mouse==grNo)&(~isnan(mean_iti_per_mouse_proficient)))...
        ,edges,bar_offset,rand_offset,'k','k',1);
    
    %                                 %Save data for  ranksum
    
    id_ii=id_ii+1;
    input_data(id_ii).data=mean_iti_per_mouse_proficient((group_no_per_mouse==grNo)&(~isnan(mean_iti_per_mouse_proficient)));
    input_data(id_ii).description=group_legend{grNo};
    
    
    
    
    bar_offset = bar_offset + 1;
    
    
    
    
    
end

title(['Mean ITI proficient per mouse'])


xticks([0 1 2])
xticklabels({'WT', 'Het', 'KO'})


ylabel('ITI (sec)')



%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for proficient ITI\n'])
[output_data] = drgMutiRanksumorTtest(input_data);

%Generate mean iti per mouse bar graph
id_ii=0;
input_data=[];

figNo = figNo +1;

try
    close(figNo)
catch
end
hFig=figure(figNo);

set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
hold on


bar_offset = 0;
edges=[0:2:50];
rand_offset=0.8;
id_ii=0;
input_data=[];

for grNo=1:max(group_no_per_mouse)
    
    switch grNo
        case 1
            bar(bar_offset,mean(mean_iti_per_mouse(group_no_per_mouse==grNo)),'g','LineWidth', 3,'EdgeColor','none')
        case 2
            bar(bar_offset,mean(mean_iti_per_mouse(group_no_per_mouse==grNo)),'b','LineWidth', 3,'EdgeColor','none')
        case 3
            bar(bar_offset,mean(mean_iti_per_mouse(group_no_per_mouse==grNo)),'y','LineWidth', 3,'EdgeColor','none')
    end
    
     
    %Violin plot
    
    [mean_out, CIout]=drgViolinPoint(mean_iti_per_mouse(group_no_per_mouse==grNo)...
        ,edges,bar_offset,rand_offset,'k','k',1);
    
    %                                 %Save data for  ranksum
    
    id_ii=id_ii+1;
    input_data(id_ii).data=mean_iti_per_mouse(group_no_per_mouse==grNo);
    input_data(id_ii).description=group_legend{grNo};
    
    
    
    
    bar_offset = bar_offset + 1;
    
    
    
    
    
end

title(['Mean ITI per mouse'])


xticks([0 1 2])
xticklabels({'WT', 'Het', 'KO'})


ylabel('ITI (sec)')



%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for  ITI\n'])
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
    
    these_mean_corr=per_corr_per_group(grNo,1:no_trials_per_group(grNo))';
    
    switch grNo
        case 1
            bar(bar_offset,mean(these_mean_corr),'b','LineWidth', 3,'EdgeColor','none')
        case 2
            bar(bar_offset,mean(these_mean_corr),'g','LineWidth', 3,'EdgeColor','none')
        case 3
            bar(bar_offset,mean(these_mean_corr),'y','LineWidth', 3,'EdgeColor','none')
    end
    
     
    %Violin plot
    
    [mean_out, CIout]=drgViolinPoint(these_mean_corr...
        ,edges,bar_offset,rand_offset,'k','k',1);
    
    %                                 %Save data for  ranksum
    
    id_ii=id_ii+1;
    input_data(id_ii).data=these_mean_corr;
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
%             [hlCR, hpCR] = boundedline([0:5:100],this_mean_per_corr_distribution', CI, 'b');
            plot([0:5:100],this_mean_per_corr_distribution', 'b')
            %             plot([0:5:100], mean_per_corr_distribution_per_group(grNo,:), '-g','LineWidth', 3)
        case 2
%             [hlCR, hpCR] = boundedline([0:5:100],this_mean_per_corr_distribution', CI, 'g');
            plot([0:5:100],this_mean_per_corr_distribution', 'g')
        case 3
%             [hlCR, hpCR] = boundedline([0:5:100],this_mean_per_corr_distribution', CI, 'y');
            plot([0:5:100],this_mean_per_corr_distribution', 'y')
    end
end
                          
%Perform the glm
fprintf(1, ['\n\nglm for percent correct distribution \n'])
tbl = table(glm_pcorr.data',glm_pcorr.per_ii',glm_pcorr.grNo',...
    'VariableNames',{'fraction','pcorr','grNo'});
mdl = fitglm(tbl,'fraction~pcorr+grNo+pcorr*grNo','CategoricalVars',[3])


id_ii=0;
input_data=[];

figNo = figNo +1;

try
    close(figNo)
catch
end
hFig=figure(figNo);

set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
hold on


bar_offset = 0;
edges=[proficient_pCorr:1:100];
rand_offset=0.8;
id_ii=0;
input_data=[];

for grNo=1:max(group_no_per_mouse)
    
    these_mean_corr1=per_corr_per_group(grNo,1:no_trials_per_group(grNo))';
    
    switch grNo
        case 1
            bar(bar_offset,mean(these_mean_corr1(these_mean_corr1>=proficient_pCorr)),'g','LineWidth', 3,'EdgeColor','none')
        case 2
            bar(bar_offset,mean(these_mean_corr1(these_mean_corr1>=proficient_pCorr)),'b','LineWidth', 3,'EdgeColor','none')
        case 3
            bar(bar_offset,mean(these_mean_corr1(these_mean_corr1>=proficient_pCorr)),'y','LineWidth', 3,'EdgeColor','none')
    end
    
     
    %Violin plot
    these_mean_corr1=these_mean_corr1(these_mean_corr1>=proficient_pCorr);
    [mean_out, CIout]=drgViolinPoint(these_mean_corr1'...
        ,edges,bar_offset,rand_offset,'k','k',1);
    
    %                                 %Save data for  ranksum
    
    id_ii=id_ii+1;
    input_data(id_ii).data=these_mean_corr1;
    input_data(id_ii).description=group_legend{grNo};
    
    
    
    
    bar_offset = bar_offset + 1;
    
    
    
    
    
end

title(['Percent correct for retreival (pCorr>' num2str(proficient_pCorr)])


xticks([0 1 2])
xticklabels({'WT', 'Het', 'KO'})


ylabel('Percent correct')

ylim([80 100])


%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for percent correct (pCorr>=' num2str(proficient_pCorr) ')\n'])
[output_data] = drgMutiRanksumorTtest(input_data);

id_ii=0;
input_data=[];

%Now plot the proficient mean percent calculated per mouse
figNo = figNo +1;

try
    close(figNo)
catch
end
hFig=figure(figNo);

set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
hold on


bar_offset = 0;
edges=[proficient_pCorr:1:100];
rand_offset=0.8;
id_ii=0;
input_data=[];

for grNo=1:max(group_no_per_mouse)
    
    these_mean_per_corr_per_mouse_prof=mean_per_corr_per_mouse_prof(group_no_per_mouse==grNo);
    these_mean_per_corr_per_mouse_prof=these_mean_per_corr_per_mouse_prof(~isnan(these_mean_per_corr_per_mouse_prof));
    switch grNo
        case 1
            bar(bar_offset,mean(these_mean_per_corr_per_mouse_prof),'g','LineWidth', 3,'EdgeColor','none')
        case 2
            bar(bar_offset,mean(these_mean_per_corr_per_mouse_prof),'b','LineWidth', 3,'EdgeColor','none')
        case 3
            bar(bar_offset,mean(these_mean_per_corr_per_mouse_prof),'y','LineWidth', 3,'EdgeColor','none')
    end
    
     
    %Violin plot
    
    [mean_out, CIout]=drgViolinPoint(these_mean_per_corr_per_mouse_prof...
        ,edges,bar_offset,rand_offset,'k','k',3);
    
    %                                 %Save data for  ranksum
    
    id_ii=id_ii+1;
    input_data(id_ii).data=these_mean_per_corr_per_mouse_prof;
    input_data(id_ii).description=group_legend{grNo};

    bar_offset = bar_offset + 1;
    

end

title(['Percent correct for retreival (pCorr>=' num2str(proficient_pCorr) ') per mouse'])


xticks([0 1 2])
xticklabels({'WT', 'Het', 'KO'})


ylabel('Percent correct')

ylim([proficient_pCorr 100])

%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for percent correct (pCorr>=' num2str(proficient_pCorr) ') per mouse\n'])
[output_data] = drgMutiRanksumorTtest(input_data);


id_ii=0;
input_data=[];

%Now plot the fraction of trials at or above proficient
figNo = figNo +1;

try
    close(figNo)
catch
end
hFig=figure(figNo);

set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
hold on


bar_offset = 0;
edges=[proficient_pCorr:1:100];
rand_offset=0.8;
id_ii=0;
input_data=[];

for grNo=1:max(group_no_per_mouse)
    
    these_fractions=100*fraction_per_corr_per_mouse_prof(group_no_per_mouse==grNo);
    these_fractions=these_fractions(~isnan(these_fractions));
    switch grNo
        case 1
            bar(bar_offset,mean(these_fractions),'g','LineWidth', 3,'EdgeColor','none')
        case 2
            bar(bar_offset,mean(these_fractions),'b','LineWidth', 3,'EdgeColor','none')
        case 3
            bar(bar_offset,mean(these_fractions),'y','LineWidth', 3,'EdgeColor','none')
    end
    
     
    %Violin plot
    
    [mean_out, CIout]=drgViolinPoint(these_fractions...
        ,edges,bar_offset,rand_offset,'k','k',3);
    
    %                                 %Save data for  ranksum
    
    id_ii=id_ii+1;
    input_data(id_ii).data=these_fractions;
    input_data(id_ii).description=group_legend{grNo};

    bar_offset = bar_offset + 1;
    

end

title(['Percent correct for retreival (pCorr>=' num2str(proficient_pCorr) ') per mouse'])


xticks([0 1 2])
xticklabels({'WT', 'Het', 'KO'})


ylabel('Fraction proficient')

ylim([20 100])

%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for percent correct (pCorr>=' num2str(proficient_pCorr) ') per mouse\n'])
[output_data] = drgMutiRanksumorTtest(input_data);

%Now plot the ratio of fraction of proficient divide by fraction of <65% 
figNo = figNo +1;

try
    close(figNo)
catch
end
hFig=figure(figNo);

set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
hold on


bar_offset = 0;
edges=[0:0.05:3];
rand_offset=0.8;
id_ii=0;
input_data=[];

for grNo=1:max(group_no_per_mouse)
    
    these_ratio_per_corr_per_mouse_prof=ratio_per_corr_per_mouse_prof(group_no_per_mouse==grNo);
    these_ratio_per_corr_per_mouse_prof=these_ratio_per_corr_per_mouse_prof((~isnan(these_ratio_per_corr_per_mouse_prof))&(~isinf(these_ratio_per_corr_per_mouse_prof)));
    switch grNo
        case 1
            bar(bar_offset,mean(these_ratio_per_corr_per_mouse_prof),'g','LineWidth', 3,'EdgeColor','none')
        case 2
            bar(bar_offset,mean(these_ratio_per_corr_per_mouse_prof),'b','LineWidth', 3,'EdgeColor','none')
        case 3
            bar(bar_offset,mean(these_ratio_per_corr_per_mouse_prof),'y','LineWidth', 3,'EdgeColor','none')
    end
    
     
    %Violin plot
    
    [mean_out, CIout]=drgViolinPoint(these_ratio_per_corr_per_mouse_prof...
        ,edges,bar_offset,rand_offset,'k','k',3);
    
    %                                 %Save data for  ranksum
    
    id_ii=id_ii+1;
    input_data(id_ii).data=these_ratio_per_corr_per_mouse_prof;
    input_data(id_ii).description=group_legend{grNo};

    bar_offset = bar_offset + 1;
    

end

title(['Ratio proficient divided by naive per mouse'])


xticks([0 1 2])
xticklabels({'WT', 'Het', 'KO'})


ylabel('Percent correct')


%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for Ratio proficient divided by naive per mouse\n'])
[output_data] = drgMutiRanksumorTtest(input_data);

%Now plot the number of sessions to proficient 
figNo = figNo +1;

try
    close(figNo)
catch
end
hFig=figure(figNo);

set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
hold on


bar_offset = 0;
edges=[0:0.05:4];
rand_offset=0.8;
id_ii=0;
input_data=[];

for grNo=1:max(group_no_per_mouse)
    
    these_sessions_to_proficient=no_sessions_to_proficient_per_mouse((group_no_per_mouse==grNo)&(no_sessions_to_proficient_per_mouse_found==1));
    
    switch grNo
        case 1
            bar(bar_offset,mean(these_sessions_to_proficient),'g','LineWidth', 3,'EdgeColor','none')
        case 2
            bar(bar_offset,mean(these_sessions_to_proficient),'b','LineWidth', 3,'EdgeColor','none')
        case 3
            bar(bar_offset,mean(these_sessions_to_proficient),'y','LineWidth', 3,'EdgeColor','none')
    end
    
     
    %Violin plot
    
    [mean_out, CIout]=drgViolinPoint(these_sessions_to_proficient...
        ,edges,bar_offset,rand_offset,'k','k',3);
    
    %                                 %Save data for  ranksum
    
    id_ii=id_ii+1;
    input_data(id_ii).data=these_sessions_to_proficient;
    input_data(id_ii).description=group_legend{grNo};

    bar_offset = bar_offset + 1;
    

end

title(['Number of sessions to proficient'])


xticks([0 1 2])
xticklabels({'WT', 'Het', 'KO'})


ylabel('Sessions')


%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for sessions to proficient\n'])
[output_data] = drgMutiRanksumorTtest(input_data);

save([matBatchPathName matFileName(1:end-4) '_' num2str(proficient_pCorr) '_out.mat'],'mean_per_corr_per_mouse_prof','group_no_per_mouse','per_corr_per_group','no_trials_per_group','ratio_per_corr_per_mouse_prof',...
    'no_sessions_to_proficient_per_mouse','no_sessions_to_proficient_per_mouse_found','mean_iti_per_mouse','mean_iti_per_mouse_proficient')
pffft=1;
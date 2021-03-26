%drgSummarizeCaMKIIbehavior
%Summarizes behavior for proficient mice
close all
clear all

group_legend{1}='WT';
group_legend{2}='Het';
group_legend{3}='KO';

%Location of files
PathName='/Users/restrepd/Documents/Projects/CaMKII_analysis/Behavior/';

%Files for 75%
FileName{1}='CaMKIICamKIAPEB_beh_03142021beh_75_out.mat';
FileName{2}='CaMKIIEAPAbeh03142021beh_75_out.mat';
FileName{3}='CaMKIIEBAP_beh_03132021beh_75_out.mat';
FileName{4}='CaMKIIPAEA_beh_03132021beh_75_out.mat';
FileName{5}='CaMKIIpz1eapa_beh_03132021beh_75_out.mat';
FileName{6}='CaMKIIpz1paea_beh_03132021beh_75_out.mat';
FileName{7}='CaMKIIpzz1ethylace_beh_03132021beh_75_out.mat';
FileName{8}='CaMKIIpzz1paea_disc_PRPall12152020beh_75_out.mat';

% %Files for 80%
% FileName{1}='CaMKIICamKIAPEB_beh_03142021beh_80_out.mat';
% FileName{2}='CaMKIIEAPAbeh03142021beh_80_out.mat';
% FileName{3}='CaMKIIEBAP_beh_03132021beh_80_out.mat';
% FileName{4}='CaMKIIPAEA_beh_03132021beh_80_out.mat';
% FileName{5}='CaMKIIpz1eapa_beh_03132021beh_80_out.mat';
% FileName{6}='CaMKIIpz1paea_beh_03132021beh_80_out.mat';
% FileName{7}='CaMKIIpzz1ethylace_beh_03132021beh_80_out.mat';
% FileName{8}='CaMKIIpzz1paea_disc_PRPall12152020beh_80_out.mat';

% %Files for 90%
% FileName{1}='CaMKIICamKIAPEB_beh_03142021beh_90_out.mat';
% FileName{2}='CaMKIIEAPAbeh03142021beh_90_out.mat';
% FileName{3}='CaMKIIEBAP_beh_03132021beh_90_out.mat';
% FileName{4}='CaMKIIPAEA_beh_03132021beh_90_out.mat';
% FileName{5}='CaMKIIpz1eapa_beh_03132021beh_90_out.mat';
% FileName{6}='CaMKIIpz1paea_beh_03132021beh_90_out.mat';
% FileName{7}='CaMKIIpzz1ethylace_beh_03132021beh_90_out.mat';
% FileName{8}='CaMKIIpzz1paea_disc_PRPall12152020beh_90_out.mat';

% %Files for 85%
% FileName{1}='CaMKIICamKIAPEB_beh_03142021beh_out.mat';
% FileName{2}='CaMKIIEAPAbeh03142021beh_out.mat';
% FileName{3}='CaMKIIEBAP_beh_03132021beh_out.mat';
% FileName{4}='CaMKIIPAEA_beh_03132021beh_out.mat';
% FileName{5}='CaMKIIpz1eapa_beh_03132021beh_out.mat';
% FileName{6}='CaMKIIpz1paea_beh_03132021beh_out.mat';
% FileName{7}='CaMKIIpzz1ethylace_beh_03132021beh_out.mat';
% FileName{8}='CaMKIIpzz1paea_disc_PRPall12152020beh_out.mat';

 
%First plot the per mouse perCorr
all_mean_per_corr_per_mouse_prof=[];
all_mean_per_odor_pair_prof=[];
all_group_no_per_mouse=[];
all_group_no_per_odor_pair=[];
ii_all=0;
ii_odor_pair=0;
glm_beh=[];
glm_ii=0; 
for filNum=1:length(FileName)
    group_no_per_mouse=[];
    mean_per_corr_per_mouse_prof=[];
    load([PathName FileName{filNum}])
    for grNo=1:3
        all_mean_per_corr_per_mouse_prof(ii_all+1:ii_all+sum(group_no_per_mouse==grNo))=mean_per_corr_per_mouse_prof(group_no_per_mouse==grNo);
        all_group_no_per_mouse(ii_all+1:ii_all+sum(group_no_per_mouse==grNo))=grNo;
        ii_all=ii_all+sum(group_no_per_mouse==grNo);
        ii_odor_pair=ii_odor_pair+1;
        all_mean_per_odor_pair_prof(ii_odor_pair)=mean(mean_per_corr_per_mouse_prof(group_no_per_mouse==grNo));
        all_group_no_per_odor_pair(ii_odor_pair)=grNo;
        
        glm_beh.data(glm_ii+1:glm_ii+sum(group_no_per_mouse==grNo))=mean_per_corr_per_mouse_prof(group_no_per_mouse==grNo);
        glm_beh.group(glm_ii+1:glm_ii+sum(group_no_per_mouse==grNo))=grNo*ones(1,sum(group_no_per_mouse==grNo));
        glm_beh.odor_pair(glm_ii+1:glm_ii+sum(group_no_per_mouse==grNo))=filNum*ones(1,sum(group_no_per_mouse==grNo));
        glm_ii=glm_ii+sum(group_no_per_mouse==grNo);
    end
end

fprintf(1, ['\n\nglm for behavior > 85% \n'])
tbl = table(glm_beh.data',glm_beh.group',glm_beh.odor_pair',...
    'VariableNames',{'pCorr','group','odor_pair'});
mdl = fitglm(tbl,'pCorr~group+odor_pair+group*odor_pair'...
    ,'CategoricalVars',[2,3])

%Now plot the proficient mean percent calculated per mouse per odor_pair
id_ii=0;
input_data=[];


figNo=0;
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
rand_offset=0.5;
id_ii=0;
input_data=[];

for grNo=1:max(all_group_no_per_mouse)
    
    these_mean_per_corr_per_mouse_prof=all_mean_per_corr_per_mouse_prof(all_group_no_per_mouse==grNo);
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

title(['Percent correct for retreival (pCorr>=80) per mouse per odor pair'])


xticks([1 2 3])
xticklabels({'WT', 'Het', 'KO'})


ylabel('Percent correct')

ylim([80 100])



%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for percent correct (pCorr>=80) per mouse per odor pair\n'])
[output_data] = drgMutiRanksumorTtest(input_data);


%Now plot the proficient mean percent calculated per odor_pair
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
edges=[80:1:100];
rand_offset=0.5;
id_ii=0;
input_data=[];

for grNo=1:max(all_group_no_per_odor_pair)
    
    these_mean_per_corr_per_mouse_prof=all_mean_per_odor_pair_prof(all_group_no_per_odor_pair==grNo);
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
    
   %Save data for  ranksum
    
    id_ii=id_ii+1;
    input_data(id_ii).data=these_mean_per_corr_per_mouse_prof;
    input_data(id_ii).description=group_legend{grNo};

    bar_offset = bar_offset + 1;
    

end

title(['Percent correct for retreival (pCorr>=80)  per odor pair'])


xticks([1 2 3])
xticklabels({'WT', 'Het', 'KO'})


ylabel('Percent correct')

ylim([80 100])



%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for percent correct (pCorr>=80) per odor pair\n'])
[output_data] = drgMutiRanksumorTtest(input_data);


pffft=1;
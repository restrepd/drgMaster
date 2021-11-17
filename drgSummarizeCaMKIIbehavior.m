%drgSummarizeCaMKIIbehavior
%Summarizes behavior for proficient mice
close all
clear all

group_legend{1}='WT';
group_legend{2}='Het';
group_legend{3}='KO';

%Location of files
PathName='/Users/restrepd/Documents/Projects/CaMKII_analysis/Behavior/';


% %Files for 85%
% pcorr_th=85;
% FileName{1}='CaMKIICamKIAPEB_beh_03142021beh_85_out.mat';
% FileName{2}='CaMKIIEAPAbeh03142021beh_85_out.mat';
% FileName{3}='CaMKIIEBAP_beh_03132021beh_85_out.mat';
% FileName{4}='CaMKIIPAEA_beh_03132021beh_85_out.mat';
% FileName{5}='CaMKIIpz1eapa_beh_03132021beh_85_out.mat';
% FileName{6}='CaMKIIpz1paea_beh_03132021beh_85_out.mat';
% FileName{7}='CaMKIIpzz1ethylace_beh_03132021beh_85_out.mat';
% FileName{8}='CaMKIIpzz1paea_beh_04022021beh_85_out.mat';

%Files for 80%
pcorr_th=80;
FileName{1}='CaMKIICamKIAPEB_beh_03142021beh_80_out.mat';
FileName{2}='CaMKIIEAPAbeh03142021beh_80_out.mat';
FileName{3}='CaMKIIEBAP_beh_03132021beh_80_out.mat';
FileName{4}='CaMKIIPAEA_beh_03132021beh_80_out.mat';
FileName{5}='CaMKIIpz1eapa_beh_03132021beh_80_out.mat';
FileName{6}='CaMKIIpz1paea_beh_03132021beh_80_out.mat';
FileName{7}='CaMKIIpzz1ethylace_beh_03132021beh_80_out.mat';
FileName{8}='CaMKIIpzz1paea_beh_04022021beh_80_out.mat';

 
%First plot the per mouse perCorr
all_mean_per_corr_per_mouse_prof=[];
all_mean_per_odor_pair_prof=[];
all_group_no_per_mouse=[];
all_group_no_per_odor_pair=[];
all_mean_iti_per_mouse=[];
all_mean_iti_per_mouse_prof=[];
all_mean_iti_per_odor_pair_prof=[];
all_mean_iti_per_odor_pair=[];
all_no_sessions_to_proficient_per_mouse=[];
ii_all=0;
ii_odor_pair=0;
glm_beh=[];
glm_ii=0;
handlesb_out.perCorr=[];
handlesb_out.mouseNo=[];
handlesb_out.odor_pairNo=[];
handlesb_out.groupNo=[];

for grNo=1:3
    for filNum=1:length(FileName)
        
        group_no_per_mouse=[];
        mean_per_corr_per_mouse_prof=[];
        
        load([PathName FileName{filNum}])
        mouseNos=[1:length(group_no_per_mouse)];
        all_mean_per_corr_per_mouse_prof(ii_all+1:ii_all+sum(group_no_per_mouse==grNo))=mean_per_corr_per_mouse_prof(group_no_per_mouse==grNo);
        all_group_no_per_mouse(ii_all+1:ii_all+sum(group_no_per_mouse==grNo))=grNo;
        all_mean_iti_per_mouse(ii_all+1:ii_all+sum(group_no_per_mouse==grNo))=mean_iti_per_mouse(group_no_per_mouse==grNo);
        all_mean_iti_per_mouse_prof(ii_all+1:ii_all+sum(group_no_per_mouse==grNo))=mean_iti_per_mouse_proficient(group_no_per_mouse==grNo);
        ii_all=ii_all+sum(group_no_per_mouse==grNo);
        
        ii_odor_pair=ii_odor_pair+1;
        all_mean_per_odor_pair_prof(ii_odor_pair)=mean(mean_per_corr_per_mouse_prof(group_no_per_mouse==grNo));
        all_mean_iti_per_odor_pair_prof(ii_odor_pair)=mean(mean_iti_per_mouse_proficient(group_no_per_mouse==grNo));
        all_mean_iti_per_odor_pair(ii_odor_pair)=mean(mean_iti_per_mouse(group_no_per_mouse==grNo));
        all_group_no_per_odor_pair(ii_odor_pair)=grNo;
        
        all_no_sessions_to_proficient_per_mouse(ii_odor_pair)=mean(no_sessions_to_proficient_per_mouse(group_no_per_mouse==grNo));
        
        glm_beh.data(glm_ii+1:glm_ii+sum(group_no_per_mouse==grNo))=mean_per_corr_per_mouse_prof(group_no_per_mouse==grNo);
        glm_beh.group(glm_ii+1:glm_ii+sum(group_no_per_mouse==grNo))=grNo*ones(1,sum(group_no_per_mouse==grNo));
        glm_beh.odor_pair(glm_ii+1:glm_ii+sum(group_no_per_mouse==grNo))=filNum*ones(1,sum(group_no_per_mouse==grNo));
        glm_ii=glm_ii+sum(group_no_per_mouse==grNo);
        
        handlesb_out.perCorr=[handlesb_out.perCorr mean_per_corr_per_mouse_prof(group_no_per_mouse==grNo)];
        handlesb_out.mouseNo=[handlesb_out.mouseNo mouseNos(group_no_per_mouse==grNo)];
        handlesb_out.odor_pairNo=[handlesb_out.odor_pairNo filNum*ones(1, sum(group_no_per_mouse==grNo))];
        handlesb_out.groupNo=[handlesb_out.groupNo grNo*ones(1, sum(group_no_per_mouse==grNo))];
    end
end

fprintf(1, ['\n\nglm for behavior > %d, per mouse per odor pair\n'],pcorr_th)
tbl = table(glm_beh.data',glm_beh.group',...
    'VariableNames',{'pCorr','group'});
mdl = fitglm(tbl,'pCorr~group'...
    ,'CategoricalVars',[2])

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

ax=gca;ax.LineWidth=3;

set(hFig, 'units','normalized','position',[.1 .5 .4 .4])
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

title(['Percent correct for retreival (pCorr>=' num2str(pcorr_th) ') per mouse per odor pair'])


xticks([0 1 2])
xticklabels({'WT', 'Het', 'KO'})


ylabel('Percent correct')

ylim([80 100])



%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for percent correct (pCorr>= %d) per mouse per odor pair\n'],pcorr_th)
[output_data] = drgMutiRanksumorTtest(input_data);


%Now plot the proficient mean over odor pairs calculated for each mouse
id_ii=0;
input_data=[];
glm_beh=[];
glm_ii=0; 

figNo = figNo +1;

try
    close(figNo)
catch
end
hFig=figure(figNo);

ax=gca;ax.LineWidth=3;

set(hFig, 'units','normalized','position',[.1 .5 .4 .4])
hold on


bar_offset = 0;
edges=[80:1:100];
rand_offset=0.5;
id_ii=0;
input_data=[];

for grNo=1:max(all_group_no_per_odor_pair)
    
    these_mean_per_corr_per_mouse_prof=all_mean_per_odor_pair_prof((all_group_no_per_odor_pair==grNo)&(~isnan(all_mean_per_odor_pair_prof)));
   
    glm_beh.data(glm_ii+1:glm_ii+length(these_mean_per_corr_per_mouse_prof))=these_mean_per_corr_per_mouse_prof;
    glm_beh.group(glm_ii+1:glm_ii+length(these_mean_per_corr_per_mouse_prof))=grNo*ones(1,length(these_mean_per_corr_per_mouse_prof));
    glm_ii=glm_ii+length(these_mean_per_corr_per_mouse_prof);
    
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

title(['Percent correct for retreival (pCorr>=' num2str(pcorr_th) ')  per mouse'])


xticks([0 1 2])
xticklabels({'WT', 'Het', 'KO'})


ylabel('Percent correct')

ylim([80 100])

fprintf(1, ['\n\nglm for behavior > %d, per mouse\n'],pcorr_th)
tbl = table(glm_beh.data',glm_beh.group',...
    'VariableNames',{'pCorr','group'});
mdl = fitglm(tbl,'pCorr~group'...
    ,'CategoricalVars',[2])


%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for percent correct (pCorr>=%d) per mouse\n'],pcorr_th)
[output_data] = drgMutiRanksumorTtest(input_data);


%Now plot the iti proficient for each mouse per odor
%pair
id_ii=0;
input_data=[];
glm_beh=[];
glm_ii=0; 

figNo = figNo +1;

try
    close(figNo)
catch
end
hFig=figure(figNo);

ax=gca;ax.LineWidth=3;

set(hFig, 'units','normalized','position',[.1 .5 .4 .4])
hold on


bar_offset = 0;
edges=[0:2:50];
rand_offset=0.5;
id_ii=0;
input_data=[];

for grNo=1:max(all_group_no_per_odor_pair)
    
    these_mean_iti_per_mouse_prof=all_mean_iti_per_mouse((all_group_no_per_mouse==grNo)&(~isnan(all_mean_iti_per_mouse)));
   
    
        
    glm_beh.data(glm_ii+1:glm_ii+length(these_mean_iti_per_mouse_prof))=these_mean_iti_per_mouse_prof;
    glm_beh.group(glm_ii+1:glm_ii+length(these_mean_iti_per_mouse_prof))=grNo*ones(1,length(these_mean_iti_per_mouse_prof));
    glm_ii=glm_ii+length(these_mean_iti_per_mouse_prof);
    
    switch grNo
        case 1
            bar(bar_offset,mean(these_mean_iti_per_mouse_prof),'g','LineWidth', 3,'EdgeColor','none')
        case 2
            bar(bar_offset,mean(these_mean_iti_per_mouse_prof),'b','LineWidth', 3,'EdgeColor','none')
        case 3
            bar(bar_offset,mean(these_mean_iti_per_mouse_prof),'y','LineWidth', 3,'EdgeColor','none')
    end
    
     
    %Violin plot
    
    [mean_out, CIout]=drgViolinPoint(these_mean_iti_per_mouse_prof...
        ,edges,bar_offset,rand_offset,'k','k',3);
    
   %Save data for  ranksum
    
    id_ii=id_ii+1;
    input_data(id_ii).data=these_mean_iti_per_mouse_prof;
    input_data(id_ii).description=group_legend{grNo};

    bar_offset = bar_offset + 1;
    

end

title(['ITI for proficient (pCorr>=' num2str(pcorr_th) ')  per mouse per odor pair'])


xticks([0 1 2])
xticklabels({'WT', 'Het', 'KO'})


ylabel('ITI')



fprintf(1, ['\n\nglm for ITI for proficient > %d, per mouse per odor pair\n'],pcorr_th)
tbl = table(glm_beh.data',glm_beh.group',...
    'VariableNames',{'pCorr','group'});
mdl = fitglm(tbl,'pCorr~group'...
    ,'CategoricalVars',[2])


%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for ITI for proficient (pCorr>=%d) per mouse per odor pair\n'],pcorr_th)
[output_data] = drgMutiRanksumorTtest(input_data);


%Now plot the iti proficient mean over odor pairs calculated for each mouse
%pair
id_ii=0;
input_data=[];
glm_beh=[];
glm_ii=0; 

figNo = figNo +1;

try
    close(figNo)
catch
end
hFig=figure(figNo);

ax=gca;ax.LineWidth=3;

set(hFig, 'units','normalized','position',[.1 .5 .4 .4])
hold on


bar_offset = 0;
edges=[0:2:50];
rand_offset=0.5;
id_ii=0;
input_data=[];

for grNo=1:max(all_group_no_per_odor_pair)
    
    these_mean_iti_per_mouse_prof=all_mean_iti_per_odor_pair((all_group_no_per_odor_pair==grNo)&(~isnan(all_mean_iti_per_odor_pair)));
   
    
        
    glm_beh.data(glm_ii+1:glm_ii+length(these_mean_iti_per_mouse_prof))=these_mean_iti_per_mouse_prof;
    glm_beh.group(glm_ii+1:glm_ii+length(these_mean_iti_per_mouse_prof))=grNo*ones(1,length(these_mean_iti_per_mouse_prof));
    glm_ii=glm_ii+length(these_mean_iti_per_mouse_prof);
    
    switch grNo
        case 1
            bar(bar_offset,mean(these_mean_iti_per_mouse_prof),'g','LineWidth', 3,'EdgeColor','none')
        case 2
            bar(bar_offset,mean(these_mean_iti_per_mouse_prof),'b','LineWidth', 3,'EdgeColor','none')
        case 3
            bar(bar_offset,mean(these_mean_iti_per_mouse_prof),'y','LineWidth', 3,'EdgeColor','none')
    end
    
     
    %Violin plot
     
    [mean_out, CIout]=drgViolinPoint(these_mean_iti_per_mouse_prof...
        ,edges,bar_offset,rand_offset,'k','k',3);
    
   %Save data for  ranksum
    
    id_ii=id_ii+1;
    input_data(id_ii).data=these_mean_iti_per_mouse_prof;
    input_data(id_ii).description=group_legend{grNo};

    bar_offset = bar_offset + 1;
    

end

title(['ITI for proficient (pCorr>=' num2str(pcorr_th) ') per mouse'])


xticks([0 1 2])
xticklabels({'WT', 'Het', 'KO'})


ylabel('ITI')



fprintf(1, ['\n\nglm for ITI for proficient > %d, per mouse\n'],pcorr_th)
tbl = table(glm_beh.data',glm_beh.group',...
    'VariableNames',{'pCorr','group'});
mdl = fitglm(tbl,'pCorr~group'...
    ,'CategoricalVars',[2])


%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for ITI for proficient (pCorr>=%d) per mouse\n'],pcorr_th)
[output_data] = drgMutiRanksumorTtest(input_data);


%Now plot the iti proficient for each mouse per odor
%pair
id_ii=0;
input_data=[];
glm_beh=[];
glm_ii=0; 

figNo = figNo +1;

try
    close(figNo)
catch
end
hFig=figure(figNo);

ax=gca;ax.LineWidth=3;

set(hFig, 'units','normalized','position',[.1 .5 .4 .4])
hold on


bar_offset = 0;
edges=[0:2:50];
rand_offset=0.5;
id_ii=0;
input_data=[];

for grNo=1:max(all_group_no_per_odor_pair)
    
    these_mean_iti_per_mouse=all_mean_iti_per_mouse((all_group_no_per_mouse==grNo)&(~isnan(all_mean_iti_per_mouse)));
   
    
        
    glm_beh.data(glm_ii+1:glm_ii+length(these_mean_iti_per_mouse))=these_mean_iti_per_mouse;
    glm_beh.group(glm_ii+1:glm_ii+length(these_mean_iti_per_mouse))=grNo*ones(1,length(these_mean_iti_per_mouse));
    glm_ii=glm_ii+length(these_mean_iti_per_mouse);
    
    switch grNo
        case 1
            bar(bar_offset,mean(these_mean_iti_per_mouse),'g','LineWidth', 3,'EdgeColor','none')
        case 2
            bar(bar_offset,mean(these_mean_iti_per_mouse),'b','LineWidth', 3,'EdgeColor','none')
        case 3
            bar(bar_offset,mean(these_mean_iti_per_mouse),'y','LineWidth', 3,'EdgeColor','none')
    end
    
     
    %Violin plot
    
    [mean_out, CIout]=drgViolinPoint(these_mean_iti_per_mouse...
        ,edges,bar_offset,rand_offset,'k','k',3);
    
   %Save data for  ranksum
    
    id_ii=id_ii+1;
    input_data(id_ii).data=these_mean_iti_per_mouse;
    input_data(id_ii).description=group_legend{grNo};

    bar_offset = bar_offset + 1;
    

end

title(['ITI for all trials (pCorr>=' num2str(pcorr_th) ')  per mouse per odor pair'])


xticks([0 1 2])
xticklabels({'WT', 'Het', 'KO'})


ylabel('ITI')



fprintf(1, ['\n\nglm for ITI for all trials > %d, per mouse per odor pair\n'],pcorr_th)
tbl = table(glm_beh.data',glm_beh.group',...
    'VariableNames',{'pCorr','group'});
mdl = fitglm(tbl,'pCorr~group'...
    ,'CategoricalVars',[2])


%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for ITI for all trials (pCorr>=%d) per mouse per odor pair\n'],pcorr_th)
[output_data] = drgMutiRanksumorTtest(input_data);


%Now plot the iti proficient mean over odor pairs calculated for each mouse
%pair
id_ii=0;
input_data=[];
glm_beh=[];
glm_ii=0; 

figNo = figNo +1;

try
    close(figNo)
catch
end
hFig=figure(figNo);

ax=gca;ax.LineWidth=3;

set(hFig, 'units','normalized','position',[.1 .5 .4 .4])
hold on


bar_offset = 0;
edges=[0:2:50];
rand_offset=0.5;
id_ii=0;
input_data=[];

for grNo=1:max(all_group_no_per_odor_pair)
    
    these_mean_iti_per_mouse=all_mean_iti_per_odor_pair((all_group_no_per_odor_pair==grNo)&(~isnan(all_mean_iti_per_odor_pair)));
   
    
        
    glm_beh.data(glm_ii+1:glm_ii+length(these_mean_iti_per_mouse))=these_mean_iti_per_mouse;
    glm_beh.group(glm_ii+1:glm_ii+length(these_mean_iti_per_mouse))=grNo*ones(1,length(these_mean_iti_per_mouse));
    glm_ii=glm_ii+length(these_mean_iti_per_mouse);
    
    switch grNo
        case 1
            bar(bar_offset,mean(these_mean_iti_per_mouse),'g','LineWidth', 3,'EdgeColor','none')
        case 2
            bar(bar_offset,mean(these_mean_iti_per_mouse),'b','LineWidth', 3,'EdgeColor','none')
        case 3
            bar(bar_offset,mean(these_mean_iti_per_mouse),'y','LineWidth', 3,'EdgeColor','none')
    end
    
     
    %Violin plot
     
    [mean_out, CIout]=drgViolinPoint(these_mean_iti_per_mouse...
        ,edges,bar_offset,rand_offset,'k','k',3);
    
   %Save data for  ranksum
    
    id_ii=id_ii+1;
    input_data(id_ii).data=these_mean_iti_per_mouse;
    input_data(id_ii).description=group_legend{grNo};

    bar_offset = bar_offset + 1;
    

end

title(['ITI for all trials (pCorr>=' num2str(pcorr_th) ') per mouse'])


xticks([0 1 2])
xticklabels({'WT', 'Het', 'KO'})


ylabel('ITI')



fprintf(1, ['\n\nglm for ITI for all trials > %d, per mouse\n'],pcorr_th)
tbl = table(glm_beh.data',glm_beh.group',...
    'VariableNames',{'pCorr','group'});
mdl = fitglm(tbl,'pCorr~group'...
    ,'CategoricalVars',[2])


%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for ITI all trials (pCorr>=%d) per mouse\n'],pcorr_th)
[output_data] = drgMutiRanksumorTtest(input_data);


%Now plot the sessions to proficiency for each mouse
id_ii=0;
input_data=[];
glm_beh=[];
glm_ii=0; 

figNo = figNo +1;

try
    close(figNo)
catch
end
hFig=figure(figNo);

ax=gca;ax.LineWidth=3;

set(hFig, 'units','normalized','position',[.1 .5 .4 .4])
hold on


bar_offset = 0;
edges=[0:2:50];
rand_offset=0.5;
id_ii=0;
input_data=[];

for grNo=1:max(all_group_no_per_odor_pair)
    
    these_mean_sessions_to_proficient=all_no_sessions_to_proficient_per_mouse((all_group_no_per_odor_pair==grNo)&(~isnan(all_mean_iti_per_odor_pair)));
   
    
        
    glm_beh.data(glm_ii+1:glm_ii+length(these_mean_sessions_to_proficient))=these_mean_sessions_to_proficient;
    glm_beh.group(glm_ii+1:glm_ii+length(these_mean_sessions_to_proficient))=grNo*ones(1,length(these_mean_sessions_to_proficient));
    glm_ii=glm_ii+length(these_mean_sessions_to_proficient);
    
    switch grNo
        case 1
            bar(bar_offset,mean(these_mean_sessions_to_proficient),'g','LineWidth', 3,'EdgeColor','none')
        case 2
            bar(bar_offset,mean(these_mean_sessions_to_proficient),'b','LineWidth', 3,'EdgeColor','none')
        case 3
            bar(bar_offset,mean(these_mean_sessions_to_proficient),'y','LineWidth', 3,'EdgeColor','none')
    end
    
     
    %Violin plot
     
    [mean_out, CIout]=drgViolinPoint(these_mean_sessions_to_proficient...
        ,edges,bar_offset,rand_offset,'k','k',5);
    
   %Save data for  ranksum
    
    id_ii=id_ii+1;
    input_data(id_ii).data=these_mean_sessions_to_proficient;
    input_data(id_ii).description=group_legend{grNo};

    bar_offset = bar_offset + 1;
    

end

title(['Sessions to proficiency (pCorr>=' num2str(pcorr_th) ')  per mouse'])


xticks([0 1 2])
xticklabels({'WT', 'Het', 'KO'})


ylabel('Sessions')



fprintf(1, ['\n\nglm for sessions to proficiency > %d, per mouse\n'],pcorr_th)
tbl = table(glm_beh.data',glm_beh.group',...
    'VariableNames',{'pCorr','group'});
mdl = fitglm(tbl,'pCorr~group'...
    ,'CategoricalVars',[2])


%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values sessions to proficiency (pCorr>=%d) per mouse\n'],pcorr_th)
[output_data] = drgMutiRanksumorTtest(input_data);

save([PathName 'drgSummarizeCaMKIIbehavior_out.mat'],'handlesb_out')

pffft=1;
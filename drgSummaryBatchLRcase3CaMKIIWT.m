function drgSummaryBatchLRcase3CaMKIIWT
%Analyzes the z scored PRP analysis performed by drgAnalysisBatchLFPCaMKIIv2 case 2
%Takes as a in input the per odorant pair .mat file output for hippocampus or prefrontal

warning('off')

close all
clear all

bandwidth_names{1}='Beta';
bandwidth_names{2}='High gamma';
bandwidth_names{3}='rSWR';

prof_naive_leg{1}='Proficient';
prof_naive_leg{2}='Naive';

group_legend{1}='WT';
group_legend{2}='Het';
group_legend{3}='KO';

evTypeLabels{1}='S+';
evTypeLabels{2}='S-';

figureNo=0;

%parameters for decision making times
decision_making_method=1;  %1 cross 0.05 followed by 0.3 sec below
%2 below critical p val followed by below for dt_crit

pCrit1=0.005;
dt_crit1=1.6666;
dt_crit2=1.2;


pCrit=log10(pCrit1);
dt_crit=0.3;
dt=0.7;
t_odor_on=0.1045;

%Location of files for proficient 80-100
PathName='/Users/restrepd/Documents/Projects/CaMKII_analysis/Lick-theta analysis/';

%Text file for statistical output
fileID = fopen([PathName 'drgSummaryBatchLRcase3CaMKIIWT.txt'],'w');

%IMPORTANT: Please note that the numbering of these files is the same as the numbering
%of odorant pairs in the mouse numbering worksheet

%Hippocampus
hippFileName{1}='spm_LFP_APEB_one_window_12102021_case3__hippocampusLFP80.mat';
hippFileName{2}='spm_LFP_EBAP_one_window_12192021_case3__hippocampusLFP80.mat';
hippFileName{3}='spm_LFP_ethyl_one_win12182021_case3__hippocampusLFP80.mat';
hippFileName{4}='spm_LFP_PAEA_one_window122121_case3__hippocampusLFP80.mat';
hippFileName{5}='spm_LFP_pz1ethyl_one_window_122321_case3__hippocampusLFP80.mat';
hippFileName{6}='spm_LFP_pz1PAEA_one_win12302021_case3__hippocampusLFP80.mat';
hippFileName{7}='spm_LFP_pzz1EAPA_one_window01012022_case3__hippocampusLFP80.mat';
hippFileName{8}='spm_LFP_pzz1PAEA_one_window12302021_case3__hippocampusLFP80.mat';
%
% %Prefrontal
% preFileName{1}='spm_LFP_APEB_one_window_12102021_case2_prefrontalLFP80.mat';
% preFileName{2}='spm_LFP_EBAP_one_window_12192021_case2_prefrontalLFP80.mat';
% preFileName{3}='spm_LFP_ethyl_one_win12182021_case2_prefrontalLFP80.mat';
% preFileName{4}='spm_LFP_PAEA_one_window122121_case2_prefrontalLFP80.mat';
% preFileName{5}='spm_LFP_pz1ethyl_one_window_122321_case2_prefrontalLFP80.mat';
% preFileName{6}='spm_LFP_pz1PAEA_one_win12302021_case2_prefrontalLFP80.mat';
% preFileName{7}='spm_LFP_pzz1EAPA_one_window01012022_case2_prefrontalLFP80.mat';
% preFileName{8}='spm_LFP_pzz1PAEA_one_window12302021_case2_prefrontalLFP80.mat';
%


%Load the table of mouse numbers
%Note: This may need to be revised for PRP
mouse_no_table='/Users/restrepd/Documents/Projects/CaMKII_analysis/Reply_to_reviewers/camkii_mice_per_odor_pair_for_PRP.xlsx';
T_mouse_no = readtable(mouse_no_table);

%Load data hippocampus
all_hippo=[];

for ii=1:length(hippFileName)
    load([PathName hippFileName{ii}])
    all_hippo(ii).handles_out=handles_out;
end


pacii=1;    %for amplitude bandwidths (beta, high gamma)


glm_licks=[];
glm_ii_licks=0;

id_ii=0;
input_data=[];

glm_licks_mm=[];
glm_ii_licks_mm=0;

id_ii_mm=0;
input_data_mm=[];

glm_from=0.5;
glm_to=2.5;


%Now plot the lick rate timecourse for naive mice
per_ii=2; %Naive

grNo=1;

%Get these z timecourses
all_SmLickRate_timecourse=[];
all_SpLickRate_timecourse=[];
these_mouse_nos=[];
ii_licks=0;
for ii=1:length(hippFileName)
    handles_out=all_hippo(ii).handles_out;
    these_mouse_no_per_op=T_mouse_no.mouse_no_per_op(T_mouse_no.odor_pair_no==ii);
    these_mouse_no=T_mouse_no.mouse_no(T_mouse_no.odor_pair_no==ii);
    
    these_mouse_numbers_for_this_op=[];
    ii_mouse=0;
    
    for lf_ii=1:handles_out.lf_ii
        if (handles_out.lickf_timecourse(lf_ii).group_no==grNo)&...
                (handles_out.lickf_timecourse(lf_ii).per_ii==per_ii)
            
            if (isfield(handles_out.lickf_timecourse(lf_ii),'Sp_lick_freq'))&(isfield(handles_out.lickf_timecourse(lf_ii),'Sm_lick_freq'))
                if (~isempty(handles_out.lickf_timecourse(lf_ii).Sp_lick_freq))&(~isempty(handles_out.lickf_timecourse(lf_ii).Sm_lick_freq))
                    ii_mouse=ii_mouse+1;
                    these_mouse_numbers_for_this_op(ii_mouse)=handles_out.lickf_timecourse(lf_ii).mouseNo;
                end
            end
        end
        
    end
    
    for ii_mouse_op=1:ii_mouse
        kk=find(these_mouse_numbers_for_this_op(ii_mouse_op)==these_mouse_no_per_op);
        these_mouse_nos(ii_mouse_op)=these_mouse_no(kk);
    end
    
    for ii_mouse=1:length(these_mouse_nos)
        ii_licks=ii_licks+1;
        
        all_SmLickRate_timecourse(ii_licks,:)=all_hippo(ii).handles_out.lick_rate.Sm_lick_freq_n(ii_mouse,:);
        all_SpLickRate_timecourse(ii_licks,:)=all_hippo(ii).handles_out.lick_rate.Sp_lick_freq_n(ii_mouse,:);
        
        all_these_mouse_nos(ii_licks)=these_mouse_nos(ii_mouse);
    end
    
end

figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.3 .3 .4 .25])


hold on

anal_t_pac=all_hippo(ii).handles_out.anal_t_pac;

CIsp = bootci(1000, @mean, all_SpLickRate_timecourse);
meansp=mean(all_SpLickRate_timecourse,1);
CIsp(1,:)=meansp-CIsp(1,:);
CIsp(2,:)=CIsp(2,:)-meansp;


CIsm = bootci(1000, @mean, all_SmLickRate_timecourse);
meansm=mean(all_SmLickRate_timecourse,1);
CIsm(1,:)=meansm-CIsm(1,:);
CIsm(2,:)=CIsm(2,:)-meansm;


if per_ii==1
    %S- Proficient
    [hlsm, hpsm] = boundedline(anal_t_pac',mean(all_SmLickRate_timecourse,1)', CIsm', 'cmap',[158/255 31/255 99/255]);
    plot(anal_t_pac',mean(all_SmLickRate_timecourse,1)','LineWidth',1,'Color',[158/255 31/255 99/255])
else
    %S- Naive
    [hlsm, hpsm] = boundedline(anal_t_pac',mean(all_SmLickRate_timecourse,1)', CIsm', 'cmap',[238/255 111/255 179/255]);
    plot(anal_t_pac',mean(all_SmLickRate_timecourse,1)','LineWidth',1,'Color',[238/255 111/255 179/255])
end

if per_ii==1
    %S+ Proficient
    [hlsp, hpsp] = boundedline(anal_t_pac',mean(all_SpLickRate_timecourse,1)', CIsp', 'cmap',[0 114/255 178/255]);
    plot(anal_t_pac',mean(all_SpLickRate_timecourse,1)','LineWidth',1,'Color',[0 114/255 178/255])
else
    %S+ naive
    [hlsp, hpsp] = boundedline(anal_t_pac',mean(all_SpLickRate_timecourse,1)', CIsp', 'cmap',[80/255 194/255 255/255]);
    plot(anal_t_pac',mean(all_SpLickRate_timecourse,1)','LineWidth',1,'Color',[80/255 194/255 255/255])
end

%S+
these_sp_data=zeros(1,size(all_SpLickRate_timecourse,1));
these_sp_data(1,:)=mean(all_SpLickRate_timecourse(:,(anal_t_pac>=glm_from)&(anal_t_pac<=glm_to)),2);
glm_licks.data(glm_ii_licks+1:glm_ii_licks+length(these_sp_data))=these_sp_data;
glm_licks.perCorr(glm_ii_licks+1:glm_ii_licks+length(these_sp_data))=per_ii*ones(1,length(these_sp_data));
glm_licks.event(glm_ii_licks+1:glm_ii_licks+length(these_sp_data))=ones(1,length(these_sp_data));
glm_ii_licks=glm_ii_licks+length(these_sp_data);

id_ii=id_ii+1;
input_data(id_ii).data=these_sp_data;
input_data(id_ii).description=['S+ ' prof_naive_leg{per_ii}];

%S-
these_sm_data=zeros(1,size(all_SmLickRate_timecourse,1));
these_sm_data(1,:)=mean(all_SmLickRate_timecourse(:,(anal_t_pac>=glm_from)&(anal_t_pac<=glm_to)),2);
glm_licks.data(glm_ii_licks+1:glm_ii_licks+length(these_sm_data))=these_sm_data;
glm_licks.perCorr(glm_ii_licks+1:glm_ii_licks+length(these_sm_data))=per_ii*ones(1,length(these_sm_data));
glm_licks.event(glm_ii_licks+1:glm_ii_licks+length(these_sm_data))=zeros(1,length(these_sm_data));
glm_ii_licks=glm_ii_licks+length(these_sm_data);

id_ii=id_ii+1;
input_data(id_ii).data=these_sm_data;
input_data(id_ii).description=['S- ' prof_naive_leg{per_ii}];

title(['Lick rate for ' prof_naive_leg{per_ii}])
ylim([0 10])
xlabel('Time(sec)')
ylabel('Hz')

%Now keep track of glm per mouse
these_unique_mouse_nos=unique(all_these_mouse_nos);
for ii_mouse=these_unique_mouse_nos
    
    %S+
    this_mm_data=mean(these_sp_data(these_unique_mouse_nos(ii_mouse)==all_these_mouse_nos));
    glm_licks_mm.data(glm_ii_licks_mm+1)=this_mm_data;
    glm_licks_mm.perCorr(glm_ii_licks_mm+1)=per_ii;
    glm_licks_mm.event(glm_ii_licks_mm+1)=1;
    glm_ii_licks_mm=glm_ii_licks_mm+1;
    
    
    input_data_mm(id_ii_mm+1).data(ii_mouse)=this_mm_data;
    
    
    %S-
    this_mm_data=mean(these_sm_data(these_unique_mouse_nos(ii_mouse)==all_these_mouse_nos));
    glm_licks_mm.data(glm_ii_licks_mm+1)=this_mm_data;
    glm_licks_mm.perCorr(glm_ii_licks_mm+1)=per_ii;
    glm_licks_mm.event(glm_ii_licks_mm+1)=0;
    glm_ii_licks_mm=glm_ii_licks_mm+1;
    
    input_data_mm(id_ii_mm+2).data=these_sm_data;
    
end
input_data_mm(id_ii_mm+1).description=['S+ ' prof_naive_leg{per_ii}];
input_data_mm(id_ii_mm+2).description=['S- ' prof_naive_leg{per_ii}];
id_ii_mm=id_ii_mm+2;

%Now plot the lick rate timecourse for proficient mice
per_ii=1; %Proficient

grNo=1;

%Get these z timecourses
all_SmLickRate_timecourse=[];
all_SpLickRate_timecourse=[];
these_mouse_nos=[];
ii_licks=0;
for ii=1:length(hippFileName)
    handles_out=all_hippo(ii).handles_out;
    these_mouse_no_per_op=T_mouse_no.mouse_no_per_op(T_mouse_no.odor_pair_no==ii);
    these_mouse_no=T_mouse_no.mouse_no(T_mouse_no.odor_pair_no==ii);
    
    these_mouse_numbers_for_this_op=[];
    ii_mouse=0;
    
    for lf_ii=1:handles_out.lf_ii
        if (handles_out.lickf_timecourse(lf_ii).group_no==grNo)&...
                (handles_out.lickf_timecourse(lf_ii).per_ii==per_ii)
            
            if (isfield(handles_out.lickf_timecourse(lf_ii),'Sp_lick_freq'))&(isfield(handles_out.lickf_timecourse(lf_ii),'Sm_lick_freq'))
                if (~isempty(handles_out.lickf_timecourse(lf_ii).Sp_lick_freq))&(~isempty(handles_out.lickf_timecourse(lf_ii).Sm_lick_freq))
                    ii_mouse=ii_mouse+1;
                    these_mouse_numbers_for_this_op(ii_mouse)=handles_out.lickf_timecourse(lf_ii).mouseNo;
                end
            end
        end
        
    end
    
    for ii_mouse_op=1:ii_mouse
        kk=find(these_mouse_numbers_for_this_op(ii_mouse_op)==these_mouse_no_per_op);
        these_mouse_nos(ii_mouse_op)=these_mouse_no(kk);
    end
    
    for ii_mouse=1:length(these_mouse_nos)
        ii_licks=ii_licks+1;
        
        all_SmLickRate_timecourse(ii_licks,:)=all_hippo(ii).handles_out.lick_rate.Sm_lick_freq_p(ii_mouse,:);
        all_SpLickRate_timecourse(ii_licks,:)=all_hippo(ii).handles_out.lick_rate.Sp_lick_freq_p(ii_mouse,:);
        
        all_these_mouse_nos(ii_licks)=these_mouse_nos(ii_mouse);
    end
    
end

figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.3 .3 .4 .25])


hold on

anal_t_pac=all_hippo(ii).handles_out.anal_t_pac;

CIsp = bootci(1000, @mean, all_SpLickRate_timecourse);
meansp=mean(all_SpLickRate_timecourse,1);
CIsp(1,:)=meansp-CIsp(1,:);
CIsp(2,:)=CIsp(2,:)-meansp;


CIsm = bootci(1000, @mean, all_SmLickRate_timecourse);
meansm=mean(all_SmLickRate_timecourse,1);
CIsm(1,:)=meansm-CIsm(1,:);
CIsm(2,:)=CIsm(2,:)-meansm;

if per_ii==1
    %S- Proficient
    [hlsm, hpsm] = boundedline(anal_t_pac',mean(all_SmLickRate_timecourse,1)', CIsm', 'cmap',[158/255 31/255 99/255]);
    plot(anal_t_pac',mean(all_SmLickRate_timecourse,1)','LineWidth',1,'Color',[158/255 31/255 99/255])
else
    %S- Naive
    [hlsm, hpsm] = boundedline(anal_t_pac',mean(all_SmLickRate_timecourse,1)', CIsm', 'cmap',[238/255 111/255 179/255]);
    plot(anal_t_pac',mean(all_SmLickRate_timecourse,1)','LineWidth',1,'Color',[238/255 111/255 179/255])
end

if per_ii==1
    %S+ Proficient
    [hlsp, hpsp] = boundedline(anal_t_pac',mean(all_SpLickRate_timecourse,1)', CIsp', 'cmap',[0 114/255 178/255]);
    plot(anal_t_pac',mean(all_SpLickRate_timecourse,1)','LineWidth',1,'Color',[0 114/255 178/255])
else
    %S+ naive
    [hlsp, hpsp] = boundedline(anal_t_pac',mean(all_SpLickRate_timecourse,1)', CIsp', 'cmap',[80/255 194/255 255/255]);
    plot(anal_t_pac',mean(all_SpLickRate_timecourse,1)','LineWidth',1,'Color',[80/255 194/255 255/255])
end

%S+
these_sp_data=zeros(1,size(all_SpLickRate_timecourse,1));
these_sp_data(1,:)=mean(all_SpLickRate_timecourse(:,(anal_t_pac>=glm_from)&(anal_t_pac<=glm_to)),2);
glm_licks.data(glm_ii_licks+1:glm_ii_licks+length(these_sp_data))=these_sp_data;
glm_licks.perCorr(glm_ii_licks+1:glm_ii_licks+length(these_sp_data))=per_ii*ones(1,length(these_sp_data));
glm_licks.event(glm_ii_licks+1:glm_ii_licks+length(these_sp_data))=ones(1,length(these_sp_data));
glm_ii_licks=glm_ii_licks+length(these_sp_data);

id_ii=id_ii+1;
input_data(id_ii).data=these_sp_data;
input_data(id_ii).description=['S+ ' prof_naive_leg{per_ii}];

%S-
these_sm_data=zeros(1,size(all_SmLickRate_timecourse,1));
these_sm_data(1,:)=mean(all_SmLickRate_timecourse(:,(anal_t_pac>=glm_from)&(anal_t_pac<=glm_to)),2);
glm_licks.data(glm_ii_licks+1:glm_ii_licks+length(these_sm_data))=these_sm_data;
glm_licks.perCorr(glm_ii_licks+1:glm_ii_licks+length(these_sm_data))=per_ii*ones(1,length(these_sm_data));
glm_licks.event(glm_ii_licks+1:glm_ii_licks+length(these_sm_data))=zeros(1,length(these_sm_data));
glm_ii_licks=glm_ii_licks+length(these_sm_data);

id_ii=id_ii+1;
input_data(id_ii).data=these_sm_data;
input_data(id_ii).description=['S- ' prof_naive_leg{per_ii}];

title(['Lick rate for ' prof_naive_leg{per_ii}])
ylim([0 10])
xlabel('Time(sec)')
ylabel('Hz')

%Now keep track of glm per mouse
these_unique_mouse_nos=unique(all_these_mouse_nos);
for ii_mouse=these_unique_mouse_nos
    
    %S+
    this_mm_data=mean(these_sp_data(these_unique_mouse_nos(ii_mouse)==all_these_mouse_nos));
    glm_licks_mm.data(glm_ii_licks_mm+1)=this_mm_data;
    glm_licks_mm.perCorr(glm_ii_licks_mm+1)=per_ii;
    glm_licks_mm.event(glm_ii_licks_mm+1)=1;
    glm_ii_licks_mm=glm_ii_licks_mm+1;
    
    
    input_data_mm(id_ii_mm+1).data(ii_mouse)=this_mm_data;
    
    
    %S-
    this_mm_data=mean(these_sm_data(these_unique_mouse_nos(ii_mouse)==all_these_mouse_nos));
    glm_licks_mm.data(glm_ii_licks_mm+1)=this_mm_data;
    glm_licks_mm.perCorr(glm_ii_licks_mm+1)=per_ii;
    glm_licks_mm.event(glm_ii_licks_mm+1)=0;
    glm_ii_licks_mm=glm_ii_licks_mm+1;
    
    input_data_mm(id_ii_mm+2).data=these_sm_data;
    
end
input_data_mm(id_ii_mm+1).description=['S+ ' prof_naive_leg{per_ii}];
input_data_mm(id_ii_mm+2).description=['S- ' prof_naive_leg{per_ii}];
id_ii_mm=id_ii_mm+2;

%Perform the glm per mouse per odor pair
fprintf(1, ['glm for lick rate per mouse per odor pair\n'])
fprintf(fileID, ['glm lick rate per mouse per odor pair\n']);

tbl = table(glm_licks.data',glm_licks.perCorr',glm_licks.event',...
    'VariableNames',{'lick_rate','naive_vs_proficient','event',});
mdl = fitglm(tbl,'lick_rate~naive_vs_proficient+event+naive_vs_proficient*event'...
    ,'CategoricalVars',[2,3])


txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);

%Do the ranksum/t-test per mouse per odor pair
fprintf(1, ['\n\nRanksum or t-test p values for lick rate  per mouse per odor pair\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values for lick rate  per mouse per odor pair\n']);

[output_data] = drgMutiRanksumorTtest(input_data, fileID);

%Perform the glm per mouse
fprintf(1, ['glm for lick rate per mouse\n'])
fprintf(fileID, ['glm lick rate per mouse\n']);

tbl = table(glm_licks_mm.data',glm_licks_mm.perCorr',glm_licks_mm.event',...
    'VariableNames',{'lick_rate','naive_vs_proficient','event',});
mdl = fitglm(tbl,'lick_rate~naive_vs_proficient+event+naive_vs_proficient*event'...
    ,'CategoricalVars',[2,3])


txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);

%Do the ranksum/t-test per mouse 
fprintf(1, ['\n\nRanksum or t-test p values for lick rate  per mouse per odor pair\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values for lick rate  per mouse per odor pair\n']);

[output_data] = drgMutiRanksumorTtest(input_data_mm, fileID);


fclose(fileID);

pffft=1;
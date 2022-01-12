function drgSummaryBatchzscorecase2CaMKIIWT
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
fileID = fopen([PathName 'drgSummaryBatchzscorecase2CaMKIIWT.txt'],'w');

%IMPORTANT: Please note that the numbering of these files is the same as the numbering
%of odorant pairs in the mouse numbering worksheet

%Hippocampus
hippFileName{1}='spm_LFP_APEB_one_window_12102021_case2__hippocampusLFP80.mat';
hippFileName{2}='spm_LFP_EBAP_one_window_12192021_case2__hippocampusLFP80.mat'; 
hippFileName{3}='spm_LFP_ethyl_one_win12182021_case2__hippocampusLFP80.mat';
hippFileName{4}='spm_LFP_PAEA_one_window122121_case2__hippocampusLFP80.mat';
hippFileName{5}='spm_LFP_pz1ethyl_one_window_122321_case2__hippocampusLFP80.mat';
hippFileName{6}='spm_LFP_pz1PAEA_one_win12302021_case2__hippocampusLFP80.mat';  
hippFileName{7}='spm_LFP_pzz1EAPA_one_window01012022_case2__hippocampusLFP80.mat'; 
hippFileName{8}='spm_LFP_pzz1PAEA_one_window12302021_case2__hippocampusLFP80.mat';

%Prefrontal
preFileName{1}='spm_LFP_APEB_one_window_12102021_case2_prefrontalLFP80.mat';
preFileName{2}='spm_LFP_EBAP_one_window_12192021_case2_prefrontalLFP80.mat'; 
preFileName{3}='spm_LFP_ethyl_one_win12182021_case2_prefrontalLFP80.mat';
preFileName{4}='spm_LFP_PAEA_one_window122121_case2_prefrontalLFP80.mat';
preFileName{5}='spm_LFP_pz1ethyl_one_window_122321_case2_prefrontalLFP80.mat';
preFileName{6}='spm_LFP_pz1PAEA_one_win12302021_case2_prefrontalLFP80.mat'; 
preFileName{7}='spm_LFP_pzz1EAPA_one_window01012022_case2_prefrontalLFP80.mat'; 
preFileName{8}='spm_LFP_pzz1PAEA_one_window12302021_case2_prefrontalLFP80.mat';



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

%Load data prefrontal
all_pre=[];

for ii=1:length(preFileName)
    load([PathName preFileName{ii}])
    all_pre(ii).handles_out=handles_out;
end

glm_dect=[];
glm_dect_ii=0;

glm_dect_mm=[];
glm_ii_mm=0;

id_dect_ii=0;
input_dect_data=[];

glm_from=0.5;
glm_to=2.5;




for pacii=[1 2]    %for amplitude bandwidths (beta, high gamma)
    
    
    glm_PRP_hipp=[];
    glm_ii_hipp=0;
    
    id_ii=0;
    input_data=[];
    
    %Now plot for the hippocampus the peak zPRP per mouse per odorant pair for S+ and S-
    for per_ii=2:-1:1
        
        grNo=1;
        
        %Get these z timecourses
        all_SmPRPtimecourse_peak=[];
        all_SpPRPtimecourse_peak=[];
        these_mouse_nos=[];
        ii_PRP=0;
        for ii=1:length(hippFileName)
            these_mouse_no_per_op=T_mouse_no.mouse_no_per_op(T_mouse_no.odor_pair_no==ii);
            these_mouse_no=T_mouse_no.mouse_no(T_mouse_no.odor_pair_no==ii);
            these_jjs=[];
            for jj=1:all_hippo(ii).handles_out.RP_ii
                if all_hippo(ii).handles_out.RPtimecourse(jj).pacii==pacii
                    if all_hippo(ii).handles_out.RPtimecourse(jj).per_ii==per_ii
                        if all_hippo(ii).handles_out.RPtimecourse(jj).group_no==grNo
                            these_jjs=[these_jjs jj];
                        end
                    end
                end
            end
            for jj=1:length(these_jjs)
                if length(all_hippo(ii).handles_out.RPtimecourse(these_jjs(jj)).meanSmPRPtimecourse_peak)>0
                    ii_PRP=ii_PRP+1;
                    all_SmPRPtimecourse_peak(ii_PRP,:)=all_hippo(ii).handles_out.RPtimecourse(these_jjs(jj)).meanSmPRPtimecourse_peak;
                    all_SpPRPtimecourse_peak(ii_PRP,:)=all_hippo(ii).handles_out.RPtimecourse(these_jjs(jj)).meanSpPRPtimecourse_peak;
                    this_nn=find(all_hippo(ii).handles_out.RPtimecourse(these_jjs(jj)).mouseNo==these_mouse_no_per_op);
                    these_mouse_nos(ii_PRP)=these_mouse_no(this_nn);
                end
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
        
        CIsp = bootci(1000, @mean, all_SpPRPtimecourse_peak);
        meansp=mean(all_SpPRPtimecourse_peak,1);
        CIsp(1,:)=meansp-CIsp(1,:);
        CIsp(2,:)=CIsp(2,:)-meansp;
        
        
        CIsm = bootci(1000, @mean, all_SmPRPtimecourse_peak);
        meansm=mean(all_SmPRPtimecourse_peak,1);
        CIsm(1,:)=meansm-CIsm(1,:);
        CIsm(2,:)=CIsm(2,:)-meansm;
        
        
        if per_ii==1
            %S- Proficient
            [hlsm, hpsm] = boundedline(anal_t_pac',mean(all_SmPRPtimecourse_peak,1)', CIsm', 'cmap',[158/255 31/255 99/255]);
        else
            %S- Naive
            [hlsm, hpsm] = boundedline(anal_t_pac',mean(all_SmPRPtimecourse_peak,1)', CIsm', 'cmap',[238/255 111/255 179/255]);
        end
        
        if per_ii==1
            %S+ Proficient
            [hlsp, hpsp] = boundedline(anal_t_pac',mean(all_SpPRPtimecourse_peak,1)', CIsp', 'cmap',[0 114/255 178/255]);
        else
            %S+ naive
            [hlsp, hpsp] = boundedline(anal_t_pac',mean(all_SpPRPtimecourse_peak,1)', CIsp', 'cmap',[80/255 194/255 255/255]);
        end
        
        %S+
        these_sp_data=zeros(1,size(all_SpPRPtimecourse_peak,1));
        these_sp_data(1,:)=mean(all_SpPRPtimecourse_peak(:,(anal_t_pac>=glm_from)&(anal_t_pac<=glm_to)),2);
        glm_PRP_hipp.data(glm_ii_hipp+1:glm_ii_hipp+length(these_sp_data))=these_sp_data;
        glm_PRP_hipp.perCorr(glm_ii_hipp+1:glm_ii_hipp+length(these_sp_data))=per_ii*ones(1,length(these_sp_data));
        glm_PRP_hipp.event(glm_ii_hipp+1:glm_ii_hipp+length(these_sp_data))=ones(1,length(these_sp_data));
        glm_PRP_hipp.peak(glm_ii_hipp+1:glm_ii_hipp+length(these_sp_data))=ones(1,length(these_sp_data));
        glm_ii_hipp=glm_ii_hipp+length(these_sp_data);
        
        id_ii=id_ii+1;
        input_data(id_ii).data=these_sp_data;
        input_data(id_ii).description=['Peak S+ ' prof_naive_leg{per_ii}];
        
        %S-
        these_sm_data=zeros(1,size(all_SmPRPtimecourse_peak,1));
        these_sm_data(1,:)=mean(all_SmPRPtimecourse_peak(:,(anal_t_pac>=glm_from)&(anal_t_pac<=glm_to)),2);
        glm_PRP_hipp.data(glm_ii_hipp+1:glm_ii_hipp+length(these_sm_data))=these_sm_data;
        glm_PRP_hipp.perCorr(glm_ii_hipp+1:glm_ii_hipp+length(these_sm_data))=per_ii*ones(1,length(these_sm_data));
        glm_PRP_hipp.event(glm_ii_hipp+1:glm_ii_hipp+length(these_sm_data))=zeros(1,length(these_sm_data));
        glm_PRP_hipp.peak(glm_ii_hipp+1:glm_ii_hipp+length(these_sm_data))=ones(1,length(these_sm_data));
        glm_ii_hipp=glm_ii_hipp+length(these_sm_data);
        
        id_ii=id_ii+1;
        input_data(id_ii).data=these_sm_data;
        input_data(id_ii).description=['Peak S- ' prof_naive_leg{per_ii}];
        
        title(['Peak zPRP for' bandwidth_names{pacii} ' ' prof_naive_leg{per_ii} ' hippocampus'])
        
        xlabel('Time(sec)')
        ylabel('zPRP')
        ylim([-2 1.3])
    end
    
    %Now plot for the hippocampus the trough zPRP per mouse per odorant pair for S+ and S-
    for per_ii=2:-1:1
        
        grNo=1;
        
        %Get these z timecourses
        all_SmPRPtimecourse_trough=[];
        all_SpPRPtimecourse_trough=[];
        ii_PRP=0;
        for ii=1:length(hippFileName)
            these_jjs=[];
            for jj=1:all_hippo(ii).handles_out.RP_ii
                if all_hippo(ii).handles_out.RPtimecourse(jj).pacii==pacii
                    if all_hippo(ii).handles_out.RPtimecourse(jj).per_ii==per_ii
                        if all_hippo(ii).handles_out.RPtimecourse(jj).group_no==grNo
                            these_jjs=[these_jjs jj];
                        end
                    end
                end
            end
            for jj=1:length(these_jjs)
                if length(all_hippo(ii).handles_out.RPtimecourse(these_jjs(jj)).meanSmPRPtimecourse_trough)>0
                    ii_PRP=ii_PRP+1;
                    all_SmPRPtimecourse_trough(ii_PRP,:)=all_hippo(ii).handles_out.RPtimecourse(these_jjs(jj)).meanSmPRPtimecourse_trough;
                    all_SpPRPtimecourse_trough(ii_PRP,:)=all_hippo(ii).handles_out.RPtimecourse(these_jjs(jj)).meanSpPRPtimecourse_trough;
                end
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
        
        CIsp = bootci(1000, @mean, all_SpPRPtimecourse_trough);
        meansp=mean(all_SpPRPtimecourse_trough,1);
        CIsp(1,:)=meansp-CIsp(1,:);
        CIsp(2,:)=CIsp(2,:)-meansp;
        
        
        CIsm = bootci(1000, @mean, all_SmPRPtimecourse_trough);
        meansm=mean(all_SmPRPtimecourse_trough,1);
        CIsm(1,:)=meansm-CIsm(1,:);
        CIsm(2,:)=CIsm(2,:)-meansm;
        
        
        if per_ii==1
            %S- Proficient
            [hlsm, hpsm] = boundedline(anal_t_pac',mean(all_SmPRPtimecourse_trough,1)', CIsm', 'cmap',[158/255 31/255 99/255]);
        else
            %S- Naive
            [hlsm, hpsm] = boundedline(anal_t_pac',mean(all_SmPRPtimecourse_trough,1)', CIsm', 'cmap',[238/255 111/255 179/255]);
        end
        
        if per_ii==1
            %S+ Proficient
            [hlsp, hpsp] = boundedline(anal_t_pac',mean(all_SpPRPtimecourse_trough,1)', CIsp', 'cmap',[0 114/255 178/255]);
        else
            %S+ naive
            [hlsp, hpsp] = boundedline(anal_t_pac',mean(all_SpPRPtimecourse_trough,1)', CIsp', 'cmap',[80/255 194/255 255/255]);
        end
        
        %S+
        these_sp_data=zeros(1,size(all_SpPRPtimecourse_trough,1));
        these_sp_data(1,:)=mean(all_SpPRPtimecourse_trough(:,(anal_t_pac>=glm_from)&(anal_t_pac<=glm_to)),2);
        glm_PRP_hipp.data(glm_ii_hipp+1:glm_ii_hipp+length(these_sp_data))=these_sp_data;
        glm_PRP_hipp.perCorr(glm_ii_hipp+1:glm_ii_hipp+length(these_sp_data))=per_ii*ones(1,length(these_sp_data));
        glm_PRP_hipp.event(glm_ii_hipp+1:glm_ii_hipp+length(these_sp_data))=ones(1,length(these_sp_data));
        glm_PRP_hipp.peak(glm_ii_hipp+1:glm_ii_hipp+length(these_sp_data))=zeros(1,length(these_sp_data));
        glm_ii_hipp=glm_ii_hipp+length(these_sp_data);
        
        id_ii=id_ii+1;
        input_data(id_ii).data=these_sp_data;
        input_data(id_ii).description=['Peak S+ ' prof_naive_leg{per_ii}];
        
        %S-
        these_sm_data=zeros(1,size(all_SmPRPtimecourse_trough,1));
        these_sm_data(1,:)=mean(all_SmPRPtimecourse_trough(:,(anal_t_pac>=glm_from)&(anal_t_pac<=glm_to)),2);
        glm_PRP_hipp.data(glm_ii_hipp+1:glm_ii_hipp+length(these_sm_data))=these_sm_data;
        glm_PRP_hipp.perCorr(glm_ii_hipp+1:glm_ii_hipp+length(these_sm_data))=per_ii*ones(1,length(these_sm_data));
        glm_PRP_hipp.event(glm_ii_hipp+1:glm_ii_hipp+length(these_sm_data))=zeros(1,length(these_sm_data));
        glm_PRP_hipp.peak(glm_ii_hipp+1:glm_ii_hipp+length(these_sm_data))=zeros(1,length(these_sm_data));
        glm_ii_hipp=glm_ii_hipp+length(these_sm_data);
        
        id_ii=id_ii+1;
        input_data(id_ii).data=these_sm_data;
        input_data(id_ii).description=['Peak S- ' prof_naive_leg{per_ii}];
        
        title(['Trough zPRP for' bandwidth_names{pacii} ' ' prof_naive_leg{per_ii} ' hippocampus'])
        
        xlabel('Time(sec)')
        ylabel('zPRP')
        ylim([-2 1.3])
    end
    
    %Perform the glm
    fprintf(1, ['glm for zPRP per mouse per odor pair for '  bandwidth_names{pacii} '\n'])
    fprintf(fileID, ['glm for zPRP per mouse per odor pair for '  bandwidth_names{pacii} '\n']);
    
    tbl = table(glm_PRP_hipp.data',glm_PRP_hipp.perCorr',glm_PRP_hipp.event',glm_PRP_hipp.peak',...
        'VariableNames',{'zPRP','naive_vs_proficient','event','peak_vs_trough'});
    mdl = fitglm(tbl,'zPRP~naive_vs_proficient+event+peak_vs_trough+naive_vs_proficient*event*peak_vs_trough'...
        ,'CategoricalVars',[2,3,4])
    
    
    txt = evalc('mdl');
    txt=regexp(txt,'<strong>','split');
    txt=cell2mat(txt);
    txt=regexp(txt,'</strong>','split');
    txt=cell2mat(txt);
    
    fprintf(fileID,'%s\n', txt);
    
    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for zPRPe per mouse per odor pair for ' bandwidth_names{pacii} ' hippocampus\n'])
    fprintf(fileID, ['\n\nRanksum or t-test p values for zPRPe per mouse per odor pair for ' bandwidth_names{pacii} ' hippocampus\n']);
    
    [output_data] = drgMutiRanksumorTtest(input_data, fileID);
    
end

%Now plot the p value timecourses for proficient mice for hippocampus
for pacii=[1 2]    %for amplitude bandwidths (beta, high gamma)
    
    grNo=1;
    
    %Get these p vals
    all_p_vals_peak=[];
    these_mouse_nos_peak=[];
    ii_peak=0;
    all_p_vals_trough=[];
    these_mouse_nos_trough=[];
    ii_trough=0;
    all_p_vals_licks=[];
    these_mouse_nos_licks=[];
    ii_licks=0;
    
    
    for ii=1:length(hippFileName)
        these_jjs=[];
        these_mouse_no_per_op=T_mouse_no.mouse_no_per_op(T_mouse_no.odor_pair_no==ii);
        these_mouse_no=T_mouse_no.mouse_no(T_mouse_no.odor_pair_no==ii);
        for jj=1:all_hippo(ii).handles_out.RP_ii
            if all_hippo(ii).handles_out.RPtimecourse(jj).pacii==pacii
                if all_hippo(ii).handles_out.RPtimecourse(jj).per_ii==per_ii
                    if all_hippo(ii).handles_out.RPtimecourse(jj).group_no==grNo
                        these_jjs=[these_jjs jj];
                    end
                end
            end
        end
        
        for jj=1:length(these_jjs)
            
            %peak
            if ~isempty(all_hippo(ii).handles_out.RPtimecourse(these_jjs(jj)).p_vals_peak)
                this_pval=all_hippo(ii).handles_out.RPtimecourse(these_jjs(jj)).p_vals_peak;
                log10_this_pval=zeros(1,length(this_pval));
                log10_this_pval(1,:)=log10(this_pval);
                log10_this_pval(isinf(log10_this_pval))=min(log10_this_pval(~isinf(log10_this_pval)));
                ii_peak=ii_peak+1;
                all_p_vals_peak(ii_peak,:)=log10_this_pval;
                this_nn=find(all_hippo(ii).handles_out.RPtimecourse(these_jjs(jj)).mouseNo==these_mouse_no_per_op);
                these_mouse_nos_peak(ii_peak)=these_mouse_no(this_nn);
            end
            
            %trough
            if ~isempty(all_hippo(ii).handles_out.RPtimecourse(these_jjs(jj)).p_vals_trough)
                this_pval=all_hippo(ii).handles_out.RPtimecourse(these_jjs(jj)).p_vals_trough;
                log10_this_pval=zeros(1,length(this_pval));
                log10_this_pval(1,:)=log10(this_pval);
                log10_this_pval(isinf(log10_this_pval))=min(log10_this_pval(~isinf(log10_this_pval)));
                ii_trough=ii_trough+1;
                all_p_vals_trough(ii_trough,:)=log10_this_pval;
                this_nn=find(all_hippo(ii).handles_out.RPtimecourse(these_jjs(jj)).mouseNo==these_mouse_no_per_op);
                these_mouse_nos_trough(ii_trough)=these_mouse_no(this_nn);
            end
            
            %licks
            if ~isempty(all_hippo(ii).handles_out.RPtimecourse(these_jjs(jj)).p_vals_licks)
                this_pval=all_hippo(ii).handles_out.RPtimecourse(these_jjs(jj)).p_vals_licks;
                if sum(isnan(this_pval))~=0
                    ii_nan=find(isnan(this_pval));
                    for ww=ii_nan
                        this_pval(ww)=this_pval(ww-1);
                    end
                end
                log10_this_pval=zeros(1,length(this_pval));
                log10_this_pval(1,:)=log10(this_pval);
                log10_this_pval(isinf(log10_this_pval))=min(log10_this_pval(~isinf(log10_this_pval)));
                ii_licks=ii_licks+1;
                all_p_vals_licks(ii_licks,:)=log10_this_pval;
                this_nn=find(all_hippo(ii).handles_out.RPtimecourse(these_jjs(jj)).mouseNo==these_mouse_no_per_op);
                these_mouse_nos_licks(ii_licks)=these_mouse_no(this_nn);
            end
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
    
    
    CIpv = bootci(1000, @mean, all_p_vals_licks);
    meanpv=mean(all_p_vals_licks,1);
    CIpv(1,:)=meanpv-CIpv(1,:);
    CIpv(2,:)=CIpv(2,:)-meanpv;
    
    
    [hlpvl, hppvl] = boundedline(anal_t_pac',mean(all_p_vals_licks,1)', CIpv', 'm');
    
    
    CIpv = bootci(1000, @mean, all_p_vals_trough);
    meanpv=mean(all_p_vals_trough,1);
    CIpv(1,:)=meanpv-CIpv(1,:);
    CIpv(2,:)=CIpv(2,:)-meanpv;
    
    
    [hlpvt, hppvt] = boundedline(anal_t_pac',mean(all_p_vals_trough,1)', CIpv', 'cmap',[213/255 94/255 0/255]);
    
    CIpv = bootci(1000, @mean, all_p_vals_peak);
    meanpv=mean(all_p_vals_peak,1);
    CIpv(1,:)=meanpv-CIpv(1,:);
    CIpv(2,:)=CIpv(2,:)-meanpv;
    
    
    [hlpvp, hppvp] = boundedline(anal_t_pac',mean(all_p_vals_peak,1)', CIpv', 'cmap',[0/255 158/255 115/255]);
    
    
    plot(anal_t_pac',mean(all_p_vals_licks,1)',  'm');
    plot(anal_t_pac',mean(all_p_vals_trough,1)',  'Color',[213/255 94/255 0/255]);
    plot(anal_t_pac',mean(all_p_vals_peak,1)',  'Color',[0/255 158/255 115/255]);
    
%     plot([anal_t_pac(1) anal_t_pac(end)],[pCrit pCrit],'-r','LineWidth', 2)
    
    ylim([-400 0])
    title(['log10 p value for ' bandwidth_names{pacii} ' for hippocampus'])
    xlabel('Time(sec)')
    ylabel('log10(p)')
    
    figureNo = figureNo + 1;
    try
        close(figureNo)
    catch
    end
    hFig=figure(figureNo);
    
    ax=gca;ax.LineWidth=3;
    set(hFig, 'units','normalized','position',[.3 .3 .3 .3])
    hold on
    
    plot(anal_t_pac',mean(all_p_vals_licks,1)',  'm');
    plot(anal_t_pac',mean(all_p_vals_trough,1)',  'Color', [213/255 94/255 0/255]);
    plot(anal_t_pac',mean(all_p_vals_peak,1)',  'Color',[0/255 158/255 115/255]);
   
    
    plot([anal_t_pac(1) anal_t_pac(end)],[pCrit pCrit],'-r','LineWidth', 2)
    xlim([-0.5 0.6])
    ylim([-40 0])
    title(['Close up of log10 p value for ' bandwidth_names{pacii} ' for hippocampus'])
    xlabel('Time(sec)')
    ylabel('log10(p)')
    
    %Find the decision making times
    dii_crit=floor(dt_crit/(anal_t_pac(2)-anal_t_pac(1)));
    dii_end=floor(dt/(anal_t_pac(2)-anal_t_pac(1)));
    
    lick_decision_times_hipp=zeros(1,size(all_p_vals_licks,1));
    found_lick_decision_times_hipp=zeros(1,size(all_p_vals_licks,1));
    lick_decision_times_hipp_control=zeros(1,size(all_p_vals_licks,1));
    found_lick_decision_times_hipp_control=zeros(1,size(all_p_vals_licks,1));
    
    %licks
    for ii_lick=1:size(all_p_vals_licks,1)
        
        this_p_val_licks=all_p_vals_licks(ii_lick,:);
        
        %Find decision time
        
        decision_making_method=1;  %1 cross 0.05 followed by 0.3 sec below
        %2 below critical p val followed by below for dt_crit
        
        %         pCrit1=0.05;
        %         dt_crit1=0.166;
        if decision_making_method==1
             %decision time
            ii_t0=find(anal_t_pac>=t_odor_on,1,'first');
            ii=ii_t0+1;
            not_found=1;
            ii_cross=ii;
            delta_ii_dt1=floor(dt_crit1/(anal_t_pac(2)-anal_t_pac(1)));
            delta_ii_dt2=ceil(dt_crit2/(anal_t_pac(2)-anal_t_pac(1)));
            
            while (ii<=ii_t0+delta_ii_dt1)&(not_found==1)
                if (this_p_val_licks(ii-1)>=log10(pCrit1))&(this_p_val_licks(ii)<=log10(pCrit1))
                    if sum(this_p_val_licks(ii:ii+delta_ii_dt2)>log10(pCrit1))==0
                        ii_cross=ii;
                        not_found=0;
                    end
                end
                ii=ii+1;
            end
            
             if ((anal_t_pac(ii_cross)+anal_t_pac(ii_cross-1))/2)-anal_t_pac(ii_t0)>dt
                not_found=1;
             end
            
            if not_found==0
                found_lick_decision_times_hipp(ii_lick)=1;
                lick_decision_times_hipp(ii_lick)=((anal_t_pac(ii_cross)+anal_t_pac(ii_cross-1))/2)-t_odor_on;
            end
            
             %control decision time
            ii_t0=find(anal_t_pac>=t_odor_on-dt,1,'first');
            ii=ii_t0+1;
            not_found=1;
            ii_cross=ii;
            delta_ii_dt1=floor(dt_crit1/(anal_t_pac(2)-anal_t_pac(1)));
            delta_ii_dt2=ceil(dt_crit2/(anal_t_pac(2)-anal_t_pac(1)));
            
            while (ii<=ii_t0+delta_ii_dt1)&(not_found==1)
                if (this_p_val_licks(ii-1)>=log10(pCrit1))&(this_p_val_licks(ii)<=log10(pCrit1))
                    if sum(this_p_val_licks(ii:ii+delta_ii_dt2)>log10(pCrit1))==0
                        ii_cross=ii;
                        not_found=0;
                    end
                end
                ii=ii+1;
            end
            
            if ((anal_t_pac(ii_cross)+anal_t_pac(ii_cross-1))/2)-anal_t_pac(ii_t0)>dt
                not_found=1;
            end
            
            if not_found==0
                found_lick_decision_times_hipp_control(ii_lick)=1;
                lick_decision_times_hipp_control(ii_lick)=((anal_t_pac(ii_cross)+anal_t_pac(ii_cross-1))/2)-t_odor_on;
            end
        else
            %decision time
            current_ii=find(anal_t_pac>=t_odor_on,1,'first');
            ii_end=current_ii+dii_end;
            
            not_found=1;
            while not_found==1
                delta_ii=find(this_p_val_licks(current_ii:end)<=pCrit,1,'first');
                if current_ii+delta_ii<=ii_end
                    if sum(this_p_val_licks(current_ii+delta_ii:current_ii+delta_ii+dii_crit)>pCrit)==0
                        not_found=0;
                        found_lick_decision_times_hipp(ii_lick)=1;
                        lick_decision_times_hipp(ii_lick)=anal_t_pac(current_ii+delta_ii)-t_odor_on;
                    else
                        current_ii=current_ii+1;
                        if current_ii>=ii_end
                            not_found=1;
                        end
                    end
                else
                    not_found=0;
                end
            end
            
            %Find control decision time
            current_ii=find(anal_t_pac>=t_odor_on-dt,1,'first');
            ii_end=current_ii+dii_end;
            
            not_found=1;
            while not_found==1
                delta_ii=find(this_p_val_licks(current_ii:end)<=pCrit,1,'first');
                if current_ii+delta_ii<=ii_end
                    if sum(this_p_val_licks(current_ii+delta_ii:current_ii+delta_ii+dii_crit)>pCrit)==0
                        not_found=0;
                        found_lick_decision_times_hipp_control(ii_lick)=1;
                        lick_decision_times_hipp_control(ii_lick)=anal_t_pac(current_ii+delta_ii)-t_odor_on;
                    else
                        current_ii=current_ii+1;
                        if current_ii>=ii_end
                            not_found=1;
                        end
                    end
                else
                    not_found=0;
                end
            end
        end
    end
    
    fprintf(1, ['For ' bandwidth_names{pacii} ' licks found %d control and %d decision times out of %d mops\n'],sum(found_lick_decision_times_hipp_control),sum(found_lick_decision_times_hipp),length(found_lick_decision_times_hipp))
    
    %peak
    peak_decision_times_hipp=zeros(1,size(all_p_vals_peak,1));
    found_peak_decision_times_hipp=zeros(1,size(all_p_vals_peak,1));
    peak_decision_times_hipp_control=zeros(1,size(all_p_vals_peak,1));
    found_peak_decision_times_hipp_control=zeros(1,size(all_p_vals_peak,1));
    
    for ii_peak=1:size(all_p_vals_peak,1)
        
        this_p_val_peak=all_p_vals_peak(ii_peak,:);
        
        if decision_making_method==1
            %decision time
            ii_t0=find(anal_t_pac>=t_odor_on,1,'first');
            ii=ii_t0+1;
            not_found=1;
            ii_cross=ii;
            delta_ii_dt1=floor(dt_crit1/(anal_t_pac(2)-anal_t_pac(1)));
            delta_ii_dt2=ceil(dt_crit2/(anal_t_pac(2)-anal_t_pac(1)));
            
            while (ii<=ii_t0+delta_ii_dt1)&(not_found==1)
                if (this_p_val_peak(ii-1)>=log10(pCrit1))&(this_p_val_peak(ii)<=log10(pCrit1))
                    if sum(this_p_val_peak(ii:ii+delta_ii_dt2)>log10(pCrit1))==0
                        ii_cross=ii;
                        not_found=0;
                    end
                end
                ii=ii+1;
            end
            
            if ((anal_t_pac(ii_cross)+anal_t_pac(ii_cross-1))/2)-anal_t_pac(ii_t0)>dt
                not_found=1;
            end
            
            if not_found==0
                found_peak_decision_times_hipp(ii_peak)=1;
                peak_decision_times_hipp(ii_peak)=((anal_t_pac(ii_cross)+anal_t_pac(ii_cross-1))/2)-t_odor_on;
            end
            
            %control decision time
            ii_t0=find(anal_t_pac>=t_odor_on-dt,1,'first');
            ii=ii_t0+1;
            not_found=1;
            ii_cross=ii;
            delta_ii_dt1=floor(dt_crit1/(anal_t_pac(2)-anal_t_pac(1)));
            delta_ii_dt2=ceil(dt_crit2/(anal_t_pac(2)-anal_t_pac(1)));
            
            while (ii<=ii_t0+delta_ii_dt1)&(not_found==1)
                if (this_p_val_peak(ii-1)>=log10(pCrit1))&(this_p_val_peak(ii)<=log10(pCrit1))
                    if sum(this_p_val_peak(ii:ii+delta_ii_dt2)>log10(pCrit1))==0
                        ii_cross=ii;
                        not_found=0;
                    end
                end
                ii=ii+1;
            end
            
            if ((anal_t_pac(ii_cross)+anal_t_pac(ii_cross-1))/2)-anal_t_pac(ii_t0)>dt
                not_found=1;
            end
            
            if not_found==0
                found_peak_decision_times_hipp_control(ii_peak)=1;
                peak_decision_times_hipp_control(ii_peak)=((anal_t_pac(ii_cross)+anal_t_pac(ii_cross-1))/2)-t_odor_on;
            end
        else
            
            %Find decision time
            current_ii=find(anal_t_pac>=t_odor_on,1,'first');
            ii_end=current_ii+dii_end;
            
            not_found=1;
            while not_found==1
                delta_ii=find(this_p_val_peak(current_ii:end)<=pCrit,1,'first');
                if current_ii+delta_ii<=ii_end
                    if sum(this_p_val_peak(current_ii+delta_ii:current_ii+delta_ii+dii_crit)>pCrit)==0
                        not_found=0;
                        found_peak_decision_times_hipp(ii_peak)=1;
                        peak_decision_times_hipp(ii_peak)=anal_t_pac(current_ii+delta_ii)-t_odor_on;
                    else
                        current_ii=current_ii+1;
                        if current_ii>=ii_end
                            not_found=1;
                        end
                    end
                else
                    not_found=0;
                end
            end
            
            %Find control decision time
            current_ii=find(anal_t_pac>=t_odor_on-dt,1,'first');
            ii_end=current_ii+dii_end;
            
            not_found=1;
            while not_found==1
                delta_ii=find(this_p_val_peak(current_ii:end)<=pCrit,1,'first');
                if current_ii+delta_ii<=ii_end
                    if sum(this_p_val_peak(current_ii+delta_ii:current_ii+delta_ii+dii_crit)>pCrit)==0
                        not_found=0;
                        found_peak_decision_times_hipp_control(ii_peak)=1;
                        peak_decision_times_hipp_control(ii_peak)=anal_t_pac(current_ii+delta_ii)-t_odor_on;
                    else
                        current_ii=current_ii+1;
                        if current_ii>=ii_end
                            not_found=1;
                        end
                    end
                else
                    not_found=0;
                end
            end
        end
    end
    
    fprintf(1, ['For ' bandwidth_names{pacii} ' peak found %d control and %d decision times out of %d mops\n'],sum(found_peak_decision_times_hipp_control),sum(found_peak_decision_times_hipp),length(found_peak_decision_times_hipp))
    
    %trough
    trough_decision_times_hipp=zeros(1,size(all_p_vals_trough,1));
    found_trough_decision_times_hipp=zeros(1,size(all_p_vals_trough,1));
    trough_decision_times_hipp_control=zeros(1,size(all_p_vals_trough,1));
    found_trough_decision_times_hipp_control=zeros(1,size(all_p_vals_trough,1));
    
    for ii_trough=1:size(all_p_vals_trough,1)
        
        this_p_val_trough=all_p_vals_trough(ii_trough,:);
        if decision_making_method==1
            %decision time
            ii_t0=find(anal_t_pac>=t_odor_on,1,'first');
            ii=ii_t0+1;
            not_found=1;
            ii_cross=ii;
            delta_ii_dt1=floor(dt_crit1/(anal_t_pac(2)-anal_t_pac(1)));
            delta_ii_dt2=ceil(dt_crit2/(anal_t_pac(2)-anal_t_pac(1)));
            
            while (ii<=ii_t0+delta_ii_dt1)&(not_found==1)
                if (this_p_val_trough(ii-1)>=log10(pCrit1))&(this_p_val_trough(ii)<=log10(pCrit1))
                    if sum(this_p_val_trough(ii:ii+delta_ii_dt2)>log10(pCrit1))==0
                        ii_cross=ii;
                        not_found=0;
                    end
                end
                ii=ii+1;
            end
            
            if ((anal_t_pac(ii_cross)+anal_t_pac(ii_cross-1))/2)-anal_t_pac(ii_t0)>dt
                not_found=1;
            end
            
            if not_found==0
                found_trough_decision_times_hipp(ii_trough)=1;
                trough_decision_times_hipp(ii_trough)=((anal_t_pac(ii_cross)+anal_t_pac(ii_cross-1))/2)-t_odor_on;
            end
            
            %control decision time
            ii_t0=find(anal_t_pac>=t_odor_on-dt,1,'first');
            ii=ii_t0+1;
            not_found=1;
            ii_cross=ii;
            delta_ii_dt1=floor(dt_crit1/(anal_t_pac(2)-anal_t_pac(1)));
            delta_ii_dt2=ceil(dt_crit2/(anal_t_pac(2)-anal_t_pac(1)));
            
            while (ii<=ii_t0+delta_ii_dt1)&(not_found==1)
                if (this_p_val_trough(ii-1)>=log10(pCrit1))&(this_p_val_trough(ii)<=log10(pCrit1))
                    if sum(this_p_val_trough(ii:ii+delta_ii_dt2)>log10(pCrit1))==0
                        ii_cross=ii;
                        not_found=0;
                    end
                end
                ii=ii+1;
            end
            
            if ((anal_t_pac(ii_cross)+anal_t_pac(ii_cross-1))/2)-anal_t_pac(ii_t0)>dt
                not_found=1;
            end
            
            if not_found==0
                found_trough_decision_times_hipp_control(ii_trough)=1;
                trough_decision_times_hipp_control(ii_trough)=((anal_t_pac(ii_cross)+anal_t_pac(ii_cross-1))/2)-t_odor_on;
            end
        else
            
            %Find decision time
            current_ii=find(anal_t_pac>=t_odor_on,1,'first');
            ii_end=current_ii+dii_end;
            
            not_found=1;
            while not_found==1
                delta_ii=find(this_p_val_trough(current_ii:end)<=pCrit,1,'first');
                if current_ii+delta_ii<=ii_end
                    if sum(this_p_val_trough(current_ii+delta_ii:current_ii+delta_ii+dii_crit)>pCrit)==0
                        not_found=0;
                        found_trough_decision_times_hipp(ii_trough)=1;
                        trough_decision_times_hipp(ii_trough)=anal_t_pac(current_ii+delta_ii)-t_odor_on;
                    else
                        current_ii=current_ii+1;
                        if current_ii>=ii_end
                            not_found=1;
                        end
                    end
                else
                    not_found=0;
                end
            end
            
            %Find control decision time
            current_ii=find(anal_t_pac>=t_odor_on-dt,1,'first');
            ii_end=current_ii+dii_end;
            
            not_found=1;
            while not_found==1
                delta_ii=find(this_p_val_trough(current_ii:end)<=pCrit,1,'first');
                if current_ii+delta_ii<=ii_end
                    if sum(this_p_val_trough(current_ii+delta_ii:current_ii+delta_ii+dii_crit)>pCrit)==0
                        not_found=0;
                        found_trough_decision_times_hipp_control(ii_trough)=1;
                        trough_decision_times_hipp_control(ii_trough)=anal_t_pac(current_ii+delta_ii)-t_odor_on;
                    else
                        current_ii=current_ii+1;
                        if current_ii>=ii_end
                            not_found=1;
                        end
                    end
                else
                    not_found=0;
                end
            end
        end
    end
    
    fprintf(1, ['For ' bandwidth_names{pacii} ' trough found %d control and %d decision times out of %d mops\n\n'],sum(found_trough_decision_times_hipp_control),sum(found_trough_decision_times_hipp),length(found_trough_decision_times_hipp))
    
    
    
    
     
    %Save data for glm for decision times
    
    %licks
    these_data=lick_decision_times_hipp(found_lick_decision_times_hipp==1);
    glm_dect.data(glm_dect_ii+1:glm_dect_ii+length(these_data))=these_data;
    glm_dect.lick_peak_trough(glm_dect_ii+1:glm_dect_ii+length(these_data))=zeros(1,length(these_data));
    glm_dect.hipp_vs_pre(glm_dect_ii+1:glm_dect_ii+length(these_data))=zeros(1,length(these_data));
    glm_dect.pacii(glm_dect_ii+1:glm_dect_ii+length(these_data))=pacii*ones(1,length(these_data));
    glm_dect_ii=glm_dect_ii+length(these_data);
    
    id_dect_ii=id_dect_ii+1;
    input_dect_data(id_dect_ii).data=these_data;
    input_dect_data(id_dect_ii).description=['Lick '];
    
    %peak
    these_data=peak_decision_times_hipp(found_peak_decision_times_hipp==1);
    glm_dect.data(glm_dect_ii+1:glm_dect_ii+length(these_data))=these_data;
    glm_dect.lick_peak_trough(glm_dect_ii+1:glm_dect_ii+length(these_data))=ones(1,length(these_data));
    glm_dect.hipp_vs_pre(glm_dect_ii+1:glm_dect_ii+length(these_data))=zeros(1,length(these_data));
    glm_dect.pacii(glm_dect_ii+1:glm_dect_ii+length(these_data))=pacii*ones(1,length(these_data));
    glm_dect_ii=glm_dect_ii+length(these_data);
    
    id_dect_ii=id_dect_ii+1;
    input_dect_data(id_dect_ii).data=these_data;
    input_dect_data(id_dect_ii).description=['Peak hippocampus ' bandwidth_names{pacii}];
    
       %trough
    these_data=trough_decision_times_hipp(found_trough_decision_times_hipp==1);
    glm_dect.data(glm_dect_ii+1:glm_dect_ii+length(these_data))=these_data;
    glm_dect.lick_peak_trough(glm_dect_ii+1:glm_dect_ii+length(these_data))=2*ones(1,length(these_data));
    glm_dect.hipp_vs_pre(glm_dect_ii+1:glm_dect_ii+length(these_data))=zeros(1,length(these_data));
    glm_dect.pacii(glm_dect_ii+1:glm_dect_ii+length(these_data))=pacii*ones(1,length(these_data));
    glm_dect_ii=glm_dect_ii+length(these_data);
    
    id_dect_ii=id_dect_ii+1;
    input_dect_data(id_dect_ii).data=these_data;
    input_dect_data(id_dect_ii).description=['Trough hippocampus ' bandwidth_names{pacii}];
    
    %Now let's plot decision times
    figureNo = figureNo + 1;
    try
        close(figureNo)
    catch
    end
    hFig=figure(figureNo);
    
    ax=gca;ax.LineWidth=3;
    set(hFig, 'units','normalized','position',[.3 .3 .3 .3])
    
    hold on
    
    all_decision_dt_peak=peak_decision_times_hipp(found_peak_decision_times_hipp==1);
    all_decision_dt_trough=trough_decision_times_hipp(found_trough_decision_times_hipp==1);
    all_decision_dt_licks=lick_decision_times_hipp(found_lick_decision_times_hipp==1);

    edges=[0:0.033:0.5];
    rand_offset=0.8;
    
    %Peak PRP
    bar_offset=1;
    bar(bar_offset,mean(all_decision_dt_peak),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
    
    
    
    %Trough PRP
    bar_offset=bar_offset+1;
    bar(bar_offset,mean(all_decision_dt_trough),'LineWidth', 3,'EdgeColor','none','FaceColor',[213/255 94/255 0/255])
    
    
    
    
    %Licks
    bar_offset=bar_offset+1;
    bar(bar_offset,mean(all_decision_dt_licks),'FaceColor','m','LineWidth', 3,'EdgeColor','none')
    
   
    %Do the mouse mean calculations and plot the mean points and lines
       
    %licks
    these_data_licks=lick_decision_times_hipp(found_lick_decision_times_hipp==1);
    these_mice_licks=these_mouse_nos_licks(found_lick_decision_times_hipp==1);
    
    %peak
    these_data_peak=peak_decision_times_hipp(found_peak_decision_times_hipp==1);
    these_mice_peak=these_mouse_nos_peak(found_peak_decision_times_hipp==1);
    
    %trough
    these_data_trough=trough_decision_times_hipp(found_trough_decision_times_hipp==1);
    these_mice_trough=these_mouse_nos_trough(found_trough_decision_times_hipp==1);
    
    unique_mouse_nos=unique([these_mice_peak these_mice_trough these_mice_licks]);
    
    for msNo=unique_mouse_nos
        
        %peak
        this_mouse_mean_peak=mean(these_data_peak(these_mice_peak==msNo));
        
        glm_dect_mm.data(glm_ii_mm+1)=this_mouse_mean_peak;
        glm_dect_mm.lick_peak_trough(glm_ii_mm+1)=0;
        glm_dect_mm.hipp_vs_pre(glm_ii_mm+1)=0;
        glm_dect_mm.pacii(glm_ii_mm+1)=pacii;
        glm_ii_mm=glm_ii_mm+1;
        
        %trough
        this_mouse_mean_trough=mean(these_data_trough(these_mice_trough==msNo));
        
        glm_dect_mm.data(glm_ii_mm+1)=this_mouse_mean_trough;
        glm_dect_mm.lick_trough_trough(glm_ii_mm+1)=1;
        glm_dect_mm.hipp_vs_pre(glm_ii_mm+1)=0;
        glm_dect_mm.pacii(glm_ii_mm+1)=pacii;
        glm_ii_mm=glm_ii_mm+1;
        
         %licks
        this_mouse_mean_licks=mean(these_data_licks(these_mice_licks==msNo));
        
        glm_dect_mm.data(glm_ii_mm+1)=this_mouse_mean_licks;
        glm_dect_mm.lick_peak_trough(glm_ii_mm+1)=2;
        glm_dect_mm.hipp_vs_pre(glm_ii_mm+1)=0;
        glm_dect_mm.pacii(glm_ii_mm+1)=pacii;
        glm_ii_mm=glm_ii_mm+1;
        
        
        plot([bar_offset-2 bar_offset-1 bar_offset],[this_mouse_mean_peak this_mouse_mean_trough this_mouse_mean_licks],'-ok','MarkerSize',6,'LineWidth',1,'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7],'Color',[0.7 0.7 0.7])
   
    end
    
    %Violin plot
    [mean_out, CIout]=drgViolinPoint(all_decision_dt_peak...
        ,edges,bar_offset-2,rand_offset,'k','k',3);
    
    %Violin plot
    [mean_out, CIout]=drgViolinPoint(all_decision_dt_trough...
        ,edges,bar_offset-1,rand_offset,'k','k',3);
    
     %Violin plot
    [mean_out, CIout]=drgViolinPoint(all_decision_dt_licks...
        ,edges,bar_offset,rand_offset,'k','k',3);
    
    
    title(['Decision time hippocampus ' bandwidth_names{pacii}])
    
    ylabel('dt')
    
    xticks([1 2 3])
    xticklabels({'Peak', 'Trough', 'Licks'})
    
    
    
end

%Now plot zPRP for prefrontal
for pacii=[1 2]    %for amplitude bandwidths (beta, high gamma)
    
    
    glm_PRP_pre=[];
    glm_ii_pre=0;
    
    id_ii=0;
    input_data=[];
    
    %Now plot for prefrontal the peak zPRP per mouse per odorant pair for S+ and S-
    for per_ii=2:-1:1
        
        grNo=1;
        
        %Get these z timecourses
        all_SmPRPtimecourse_peak=[];
        all_SpPRPtimecourse_peak=[];
        these_mouse_nos=[];
        ii_PRP=0;
        for ii=1:length(preFileName)
            these_mouse_no_per_op=T_mouse_no.mouse_no_per_op(T_mouse_no.odor_pair_no==ii);
            these_mouse_no=T_mouse_no.mouse_no(T_mouse_no.odor_pair_no==ii);
            these_jjs=[];
            for jj=1:all_pre(ii).handles_out.RP_ii
                if all_pre(ii).handles_out.RPtimecourse(jj).pacii==pacii
                    if all_pre(ii).handles_out.RPtimecourse(jj).per_ii==per_ii
                        if all_pre(ii).handles_out.RPtimecourse(jj).group_no==grNo
                            these_jjs=[these_jjs jj];
                        end
                    end
                end
            end
            for jj=1:length(these_jjs)
                if length(all_pre(ii).handles_out.RPtimecourse(these_jjs(jj)).meanSmPRPtimecourse_peak)>0
                    ii_PRP=ii_PRP+1;
                    all_SmPRPtimecourse_peak(ii_PRP,:)=all_pre(ii).handles_out.RPtimecourse(these_jjs(jj)).meanSmPRPtimecourse_peak;
                    all_SpPRPtimecourse_peak(ii_PRP,:)=all_pre(ii).handles_out.RPtimecourse(these_jjs(jj)).meanSpPRPtimecourse_peak;
                    this_nn=find(all_pre(ii).handles_out.RPtimecourse(these_jjs(jj)).mouseNo==these_mouse_no_per_op);
                    these_mouse_nos(ii_PRP)=these_mouse_no(this_nn);
                end
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
        
        anal_t_pac=all_pre(ii).handles_out.anal_t_pac;
        
        CIsp = bootci(1000, @mean, all_SpPRPtimecourse_peak);
        meansp=mean(all_SpPRPtimecourse_peak,1);
        CIsp(1,:)=meansp-CIsp(1,:);
        CIsp(2,:)=CIsp(2,:)-meansp;
        
        
        CIsm = bootci(1000, @mean, all_SmPRPtimecourse_peak);
        meansm=mean(all_SmPRPtimecourse_peak,1);
        CIsm(1,:)=meansm-CIsm(1,:);
        CIsm(2,:)=CIsm(2,:)-meansm;
        
        
        if per_ii==1
            %S- Proficient
            [hlsm, hpsm] = boundedline(anal_t_pac',mean(all_SmPRPtimecourse_peak,1)', CIsm', 'cmap',[158/255 31/255 99/255]);
        else
            %S- Naive
            [hlsm, hpsm] = boundedline(anal_t_pac',mean(all_SmPRPtimecourse_peak,1)', CIsm', 'cmap',[238/255 111/255 179/255]);
        end
        
        if per_ii==1
            %S+ Proficient
            [hlsp, hpsp] = boundedline(anal_t_pac',mean(all_SpPRPtimecourse_peak,1)', CIsp', 'cmap',[0 114/255 178/255]);
        else
            %S+ naive
            [hlsp, hpsp] = boundedline(anal_t_pac',mean(all_SpPRPtimecourse_peak,1)', CIsp', 'cmap',[80/255 194/255 255/255]);
        end
        
        %S+
        these_sp_data=zeros(1,size(all_SpPRPtimecourse_peak,1));
        these_sp_data(1,:)=mean(all_SpPRPtimecourse_peak(:,(anal_t_pac>=glm_from)&(anal_t_pac<=glm_to)),2);
        glm_PRP_pre.data(glm_ii_pre+1:glm_ii_pre+length(these_sp_data))=these_sp_data;
        glm_PRP_pre.perCorr(glm_ii_pre+1:glm_ii_pre+length(these_sp_data))=per_ii*ones(1,length(these_sp_data));
        glm_PRP_pre.event(glm_ii_pre+1:glm_ii_pre+length(these_sp_data))=ones(1,length(these_sp_data));
        glm_PRP_pre.peak(glm_ii_pre+1:glm_ii_pre+length(these_sp_data))=ones(1,length(these_sp_data));
        glm_ii_pre=glm_ii_pre+length(these_sp_data);
        
        id_ii=id_ii+1;
        input_data(id_ii).data=these_sp_data;
        input_data(id_ii).description=['Peak S+ ' prof_naive_leg{per_ii}];
        
        %S-
        these_sm_data=zeros(1,size(all_SmPRPtimecourse_peak,1));
        these_sm_data(1,:)=mean(all_SmPRPtimecourse_peak(:,(anal_t_pac>=glm_from)&(anal_t_pac<=glm_to)),2);
        glm_PRP_pre.data(glm_ii_pre+1:glm_ii_pre+length(these_sm_data))=these_sm_data;
        glm_PRP_pre.perCorr(glm_ii_pre+1:glm_ii_pre+length(these_sm_data))=per_ii*ones(1,length(these_sm_data));
        glm_PRP_pre.event(glm_ii_pre+1:glm_ii_pre+length(these_sm_data))=zeros(1,length(these_sm_data));
        glm_PRP_pre.peak(glm_ii_pre+1:glm_ii_pre+length(these_sm_data))=ones(1,length(these_sm_data));
        glm_ii_pre=glm_ii_pre+length(these_sm_data);
        
        id_ii=id_ii+1;
        input_data(id_ii).data=these_sm_data;
        input_data(id_ii).description=['Peak S- ' prof_naive_leg{per_ii}];
        
        title(['Peak zPRP for' bandwidth_names{pacii} ' ' prof_naive_leg{per_ii} ' prefrontal'])
        
        xlabel('Time(sec)')
        ylabel('zPRP')
        ylim([-2 1.3])
    end
    
    %Now plot for prefrontal the trough zPRP per mouse per odorant pair for S+ and S-
    for per_ii=2:-1:1
        
        grNo=1;
        
        %Get these z timecourses
        all_SmPRPtimecourse_trough=[];
        all_SpPRPtimecourse_trough=[];
        ii_PRP=0;
        for ii=1:length(preFileName)
            these_jjs=[];
            for jj=1:all_pre(ii).handles_out.RP_ii
                if all_pre(ii).handles_out.RPtimecourse(jj).pacii==pacii
                    if all_pre(ii).handles_out.RPtimecourse(jj).per_ii==per_ii
                        if all_pre(ii).handles_out.RPtimecourse(jj).group_no==grNo
                            these_jjs=[these_jjs jj];
                        end
                    end
                end
            end
            for jj=1:length(these_jjs)
                if length(all_pre(ii).handles_out.RPtimecourse(these_jjs(jj)).meanSmPRPtimecourse_trough)>0
                    ii_PRP=ii_PRP+1;
                    all_SmPRPtimecourse_trough(ii_PRP,:)=all_pre(ii).handles_out.RPtimecourse(these_jjs(jj)).meanSmPRPtimecourse_trough;
                    all_SpPRPtimecourse_trough(ii_PRP,:)=all_pre(ii).handles_out.RPtimecourse(these_jjs(jj)).meanSpPRPtimecourse_trough;
                end
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
        
        anal_t_pac=all_pre(ii).handles_out.anal_t_pac;
        
        CIsp = bootci(1000, @mean, all_SpPRPtimecourse_trough);
        meansp=mean(all_SpPRPtimecourse_trough,1);
        CIsp(1,:)=meansp-CIsp(1,:);
        CIsp(2,:)=CIsp(2,:)-meansp;
        
        
        CIsm = bootci(1000, @mean, all_SmPRPtimecourse_trough);
        meansm=mean(all_SmPRPtimecourse_trough,1);
        CIsm(1,:)=meansm-CIsm(1,:);
        CIsm(2,:)=CIsm(2,:)-meansm;
        
        
        if per_ii==1
            %S- Proficient
            [hlsm, hpsm] = boundedline(anal_t_pac',mean(all_SmPRPtimecourse_trough,1)', CIsm', 'cmap',[158/255 31/255 99/255]);
        else
            %S- Naive
            [hlsm, hpsm] = boundedline(anal_t_pac',mean(all_SmPRPtimecourse_trough,1)', CIsm', 'cmap',[238/255 111/255 179/255]);
        end
        
        if per_ii==1
            %S+ Proficient
            [hlsp, hpsp] = boundedline(anal_t_pac',mean(all_SpPRPtimecourse_trough,1)', CIsp', 'cmap',[0 114/255 178/255]);
        else
            %S+ naive
            [hlsp, hpsp] = boundedline(anal_t_pac',mean(all_SpPRPtimecourse_trough,1)', CIsp', 'cmap',[80/255 194/255 255/255]);
        end
        
        %S+
        these_sp_data=zeros(1,size(all_SpPRPtimecourse_trough,1));
        these_sp_data(1,:)=mean(all_SpPRPtimecourse_trough(:,(anal_t_pac>=glm_from)&(anal_t_pac<=glm_to)),2);
        glm_PRP_pre.data(glm_ii_pre+1:glm_ii_pre+length(these_sp_data))=these_sp_data;
        glm_PRP_pre.perCorr(glm_ii_pre+1:glm_ii_pre+length(these_sp_data))=per_ii*ones(1,length(these_sp_data));
        glm_PRP_pre.event(glm_ii_pre+1:glm_ii_pre+length(these_sp_data))=ones(1,length(these_sp_data));
        glm_PRP_pre.peak(glm_ii_pre+1:glm_ii_pre+length(these_sp_data))=zeros(1,length(these_sp_data));
        glm_ii_pre=glm_ii_pre+length(these_sp_data);
        
        id_ii=id_ii+1;
        input_data(id_ii).data=these_sp_data;
        input_data(id_ii).description=['Peak S+ ' prof_naive_leg{per_ii}];
        
        %S-
        these_sm_data=zeros(1,size(all_SmPRPtimecourse_trough,1));
        these_sm_data(1,:)=mean(all_SmPRPtimecourse_trough(:,(anal_t_pac>=glm_from)&(anal_t_pac<=glm_to)),2);
        glm_PRP_pre.data(glm_ii_pre+1:glm_ii_pre+length(these_sm_data))=these_sm_data;
        glm_PRP_pre.perCorr(glm_ii_pre+1:glm_ii_pre+length(these_sm_data))=per_ii*ones(1,length(these_sm_data));
        glm_PRP_pre.event(glm_ii_pre+1:glm_ii_pre+length(these_sm_data))=zeros(1,length(these_sm_data));
        glm_PRP_pre.peak(glm_ii_pre+1:glm_ii_pre+length(these_sm_data))=zeros(1,length(these_sm_data));
        glm_ii_pre=glm_ii_pre+length(these_sm_data);
        
        id_ii=id_ii+1;
        input_data(id_ii).data=these_sm_data;
        input_data(id_ii).description=['Peak S- ' prof_naive_leg{per_ii}];
        
        title(['Trough zPRP for' bandwidth_names{pacii} ' ' prof_naive_leg{per_ii} ' prefrontal'])
        
        xlabel('Time(sec)')
        ylabel('zPRP')
        ylim([-2 1.3])
    end
    
    %Perform the glm
    fprintf(1, ['glm for zPRP per mouse per odor pair for '  bandwidth_names{pacii} ' prefrontal\n'])
    fprintf(fileID, ['glm for zPRP per mouse per odor pair for '  bandwidth_names{pacii} ' prefrontal\n']);
    
    tbl = table(glm_PRP_pre.data',glm_PRP_pre.perCorr',glm_PRP_pre.event',glm_PRP_pre.peak',...
        'VariableNames',{'zPRP','naive_vs_proficient','event','peak_vs_trough'});
    mdl = fitglm(tbl,'zPRP~naive_vs_proficient+event+peak_vs_trough+naive_vs_proficient*event*peak_vs_trough'...
        ,'CategoricalVars',[2,3,4])
    
    
    txt = evalc('mdl');
    txt=regexp(txt,'<strong>','split');
    txt=cell2mat(txt);
    txt=regexp(txt,'</strong>','split');
    txt=cell2mat(txt);
    
    fprintf(fileID,'%s\n', txt);
    
    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for zPRPe per mouse per odor pair for ' bandwidth_names{pacii} ' prefrontal\n'])
    fprintf(fileID, ['\n\nRanksum or t-test p values for zPRPe per mouse per odor pair for ' bandwidth_names{pacii} ' prefrontal\n']);
    
    [output_data] = drgMutiRanksumorTtest(input_data, fileID);
    
end

%Now plot the p value timecourses for proficient mice for prefrontal
for pacii=[1 2]    %for amplitude bandwidths (beta, high gamma)
    
    grNo=1;
    
    %Get these p vals
    all_p_vals_peak=[];
    ii_peak=0;
    all_p_vals_trough=[];
    ii_trough=0;
    all_p_vals_licks=[];
    ii_licks=0;
    
    for ii=1:length(preFileName)
        these_jjs=[];
        for jj=1:all_pre(ii).handles_out.RP_ii
            if all_pre(ii).handles_out.RPtimecourse(jj).pacii==pacii
                if all_pre(ii).handles_out.RPtimecourse(jj).per_ii==per_ii
                    if all_pre(ii).handles_out.RPtimecourse(jj).group_no==grNo
                        these_jjs=[these_jjs jj];
                    end
                end
            end
        end
        
        for jj=1:length(these_jjs)
            
            %peak
            if ~isempty(all_pre(ii).handles_out.RPtimecourse(these_jjs(jj)).p_vals_peak)
                this_pval=all_pre(ii).handles_out.RPtimecourse(these_jjs(jj)).p_vals_peak;
                log10_this_pval=zeros(1,length(this_pval));
                log10_this_pval(1,:)=log10(this_pval);
                log10_this_pval(isinf(log10_this_pval))=min(log10_this_pval(~isinf(log10_this_pval)));
                ii_peak=ii_peak+1;
                all_p_vals_peak(ii_peak,:)=log10_this_pval;
            end
            
            %trough
            if ~isempty(all_pre(ii).handles_out.RPtimecourse(these_jjs(jj)).p_vals_trough)
                this_pval=all_pre(ii).handles_out.RPtimecourse(these_jjs(jj)).p_vals_trough;
                log10_this_pval=zeros(1,length(this_pval));
                log10_this_pval(1,:)=log10(this_pval);
                log10_this_pval(isinf(log10_this_pval))=min(log10_this_pval(~isinf(log10_this_pval)));
                ii_trough=ii_trough+1;
                all_p_vals_trough(ii_trough,:)=log10_this_pval;
            end
            
            %licks
            if ~isempty(all_pre(ii).handles_out.RPtimecourse(these_jjs(jj)).p_vals_licks)
                this_pval=all_pre(ii).handles_out.RPtimecourse(these_jjs(jj)).p_vals_licks;
                if sum(isnan(this_pval))~=0
                    ii_nan=find(isnan(this_pval));
                    for jj=ii_nan
                        this_pval(jj)=this_pval(jj-1);
                    end
                end
                log10_this_pval=zeros(1,length(this_pval));
                log10_this_pval(1,:)=log10(this_pval);
                log10_this_pval(isinf(log10_this_pval))=min(log10_this_pval(~isinf(log10_this_pval)));
                ii_licks=ii_licks+1;
                all_p_vals_licks(ii_licks,:)=log10_this_pval;
            end
            
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
    

    CIpv = bootci(1000, @mean, all_p_vals_licks);
    meanpv=mean(all_p_vals_licks,1);
    CIpv(1,:)=meanpv-CIpv(1,:);
    CIpv(2,:)=CIpv(2,:)-meanpv;
    
    
    [hlpvl, hppvl] = boundedline(anal_t_pac',mean(all_p_vals_licks,1)', CIpv', 'm');
    
    
    CIpv = bootci(1000, @mean, all_p_vals_trough);
    meanpv=mean(all_p_vals_trough,1);
    CIpv(1,:)=meanpv-CIpv(1,:);
    CIpv(2,:)=CIpv(2,:)-meanpv;
    
    
    [hlpvt, hppvt] = boundedline(anal_t_pac',mean(all_p_vals_trough,1)', CIpv', 'cmap',[213/255 94/255 0/255]);
    
    CIpv = bootci(1000, @mean, all_p_vals_peak);
    meanpv=mean(all_p_vals_peak,1);
    CIpv(1,:)=meanpv-CIpv(1,:);
    CIpv(2,:)=CIpv(2,:)-meanpv;
    
    
    [hlpvp, hppvp] = boundedline(anal_t_pac',mean(all_p_vals_peak,1)', CIpv', 'cmap',[0/255 158/255 115/255]);
    
    
    plot(anal_t_pac',mean(all_p_vals_licks,1)',  'm');
    plot(anal_t_pac',mean(all_p_vals_trough,1)',  'Color', [213/255 94/255 0/255]);
    plot(anal_t_pac',mean(all_p_vals_peak,1)',  'Color',[0/255 158/255 115/255]);
    
%     plot([anal_t_pac(1) anal_t_pac(end)],[pCrit pCrit],'-r','LineWidth', 2)
    
    ylim([-400 0])
    title(['log10 p value for ' bandwidth_names{pacii} ' for prefrontal'])
    xlabel('Time(sec)')
    ylabel('log10(p)')
    
    figureNo = figureNo + 1;
    try
        close(figureNo)
    catch
    end
    hFig=figure(figureNo);
    
    ax=gca;ax.LineWidth=3;
    set(hFig, 'units','normalized','position',[.3 .3 .3 .3])
    hold on
    
    plot(anal_t_pac',mean(all_p_vals_licks,1)',  'm');
    plot(anal_t_pac',mean(all_p_vals_trough,1)',  'Color', [213/255 94/255 0/255]);
    plot(anal_t_pac',mean(all_p_vals_peak,1)',  'Color',[0/255 158/255 115/255]);
    
  
    
    plot([anal_t_pac(1) anal_t_pac(end)],[pCrit pCrit],'-r','LineWidth', 2)
    xlim([-0.5 0.6])
    ylim([-40 0])
    title(['Close up of log10 p value for ' bandwidth_names{pacii} ' for prefrontal'])
    xlabel('Time(sec)')
    ylabel('log10(p)')
    %Find the decision making times
    dii_crit=floor(dt_crit/(anal_t_pac(2)-anal_t_pac(1)));
    dii_end=floor(dt/(anal_t_pac(2)-anal_t_pac(1)));
    

    %peak
    peak_decision_times_pre=zeros(1,size(all_p_vals_peak,1));
    found_peak_decision_times_pre=zeros(1,size(all_p_vals_peak,1));
    peak_decision_times_pre_control=zeros(1,size(all_p_vals_peak,1));
    found_peak_decision_times_pre_control=zeros(1,size(all_p_vals_peak,1));
    
    for ii_peak=1:size(all_p_vals_peak,1)
        
        this_p_val_peak=all_p_vals_peak(ii_peak,:);
        
        if decision_making_method==1
            %decision time
            ii_t0=find(anal_t_pac>=t_odor_on,1,'first');
            ii=ii_t0+1;
            not_found=1;
            ii_cross=ii;
            delta_ii_dt1=floor(dt_crit1/(anal_t_pac(2)-anal_t_pac(1)));
            delta_ii_dt2=ceil(dt_crit2/(anal_t_pac(2)-anal_t_pac(1)));
            
            while (ii<=ii_t0+delta_ii_dt1)&(not_found==1)
                if (this_p_val_peak(ii-1)>=log10(pCrit1))&(this_p_val_peak(ii)<=log10(pCrit1))
                    if sum(this_p_val_peak(ii:ii+delta_ii_dt2)>log10(pCrit1))==0
                        ii_cross=ii;
                        not_found=0;
                    end
                end
                ii=ii+1;
            end
            
            if ((anal_t_pac(ii_cross)+anal_t_pac(ii_cross-1))/2)-anal_t_pac(ii_t0)>dt
                not_found=1;
            end
            
            if not_found==0
                found_peak_decision_times_pre(ii_peak)=1;
                peak_decision_times_pre(ii_peak)=((anal_t_pac(ii_cross)+anal_t_pac(ii_cross-1))/2)-t_odor_on;
            end
            
            %control decision time
            ii_t0=find(anal_t_pac>=t_odor_on-dt,1,'first');
            ii=ii_t0+1;
            not_found=1;
            ii_cross=ii;
            delta_ii_dt1=floor(dt_crit1/(anal_t_pac(2)-anal_t_pac(1)));
            delta_ii_dt2=ceil(dt_crit2/(anal_t_pac(2)-anal_t_pac(1)));
            
            while (ii<=ii_t0+delta_ii_dt1)&(not_found==1)
                if (this_p_val_peak(ii-1)>=log10(pCrit1))&(this_p_val_peak(ii)<=log10(pCrit1))
                    if sum(this_p_val_peak(ii:ii+delta_ii_dt2)>log10(pCrit1))==0
                        ii_cross=ii;
                        not_found=0;
                    end
                end
                ii=ii+1;
            end
            
            if ((anal_t_pac(ii_cross)+anal_t_pac(ii_cross-1))/2)-anal_t_pac(ii_t0)>dt
                not_found=1;
            end
            
            if not_found==0
                found_peak_decision_times_pre_control(ii_peak)=1;
                peak_decision_times_pre_control(ii_peak)=((anal_t_pac(ii_cross)+anal_t_pac(ii_cross-1))/2)-t_odor_on;
            end
        else
            
            %Find decision time
            current_ii=find(anal_t_pac>=t_odor_on,1,'first');
            ii_end=current_ii+dii_end;
            
            not_found=1;
            while not_found==1
                delta_ii=find(this_p_val_peak(current_ii:end)<=pCrit,1,'first');
                if current_ii+delta_ii<=ii_end
                    if sum(this_p_val_peak(current_ii+delta_ii:current_ii+delta_ii+dii_crit)>pCrit)==0
                        not_found=0;
                        found_peak_decision_times_pre(ii_peak)=1;
                        peak_decision_times_pre(ii_peak)=anal_t_pac(current_ii+delta_ii)-t_odor_on;
                    else
                        current_ii=current_ii+1;
                        if current_ii>=ii_end
                            not_found=1;
                        end
                    end
                else
                    not_found=0;
                end
            end
            
            %Find control decision time
            current_ii=find(anal_t_pac>=t_odor_on-dt,1,'first');
            ii_end=current_ii+dii_end;
            
            not_found=1;
            while not_found==1
                delta_ii=find(this_p_val_peak(current_ii:end)<=pCrit,1,'first');
                if current_ii+delta_ii<=ii_end
                    if sum(this_p_val_peak(current_ii+delta_ii:current_ii+delta_ii+dii_crit)>pCrit)==0
                        not_found=0;
                        found_peak_decision_times_pre_control(ii_peak)=1;
                        peak_decision_times_pre_control(ii_peak)=anal_t_pac(current_ii+delta_ii)-t_odor_on;
                    else
                        current_ii=current_ii+1;
                        if current_ii>=ii_end
                            not_found=1;
                        end
                    end
                else
                    not_found=0;
                end
            end
        end
    end
    
    fprintf(1, ['For ' bandwidth_names{pacii} ' peak found %d control and %d decision times out of %d mops\n'],sum(found_peak_decision_times_pre_control),sum(found_peak_decision_times_pre),length(found_peak_decision_times_pre))
    
    %trough
    trough_decision_times_pre=zeros(1,size(all_p_vals_trough,1));
    found_trough_decision_times_pre=zeros(1,size(all_p_vals_trough,1));
    trough_decision_times_pre_control=zeros(1,size(all_p_vals_trough,1));
    found_trough_decision_times_pre_control=zeros(1,size(all_p_vals_trough,1));
    
    for ii_trough=1:size(all_p_vals_trough,1)
        
        this_p_val_trough=all_p_vals_trough(ii_trough,:);
        if decision_making_method==1
            %decision time
            ii_t0=find(anal_t_pac>=t_odor_on,1,'first');
            ii=ii_t0+1;
            not_found=1;
            ii_cross=ii;
            delta_ii_dt1=floor(dt_crit1/(anal_t_pac(2)-anal_t_pac(1)));
            delta_ii_dt2=ceil(dt_crit2/(anal_t_pac(2)-anal_t_pac(1)));
            
            while (ii<=ii_t0+delta_ii_dt1)&(not_found==1)
                if (this_p_val_trough(ii-1)>=log10(pCrit1))&(this_p_val_trough(ii)<=log10(pCrit1))
                    if sum(this_p_val_trough(ii:ii+delta_ii_dt2)>log10(pCrit1))==0
                        ii_cross=ii;
                        not_found=0;
                    end
                end
                ii=ii+1;
            end
            
            if ((anal_t_pac(ii_cross)+anal_t_pac(ii_cross-1))/2)-anal_t_pac(ii_t0)>dt
                not_found=1;
            end
            
            if not_found==0
                found_trough_decision_times_pre(ii_trough)=1;
                trough_decision_times_pre(ii_trough)=((anal_t_pac(ii_cross)+anal_t_pac(ii_cross-1))/2)-t_odor_on;
            end
            
            %control decision time
            ii_t0=find(anal_t_pac>=t_odor_on-dt,1,'first');
            ii=ii_t0+1;
            not_found=1;
            ii_cross=ii;
            delta_ii_dt1=floor(dt_crit1/(anal_t_pac(2)-anal_t_pac(1)));
            delta_ii_dt2=ceil(dt_crit2/(anal_t_pac(2)-anal_t_pac(1)));
            
            while (ii<=ii_t0+delta_ii_dt1)&(not_found==1)
                if (this_p_val_trough(ii-1)>=log10(pCrit1))&(this_p_val_trough(ii)<=log10(pCrit1))
                    if sum(this_p_val_trough(ii:ii+delta_ii_dt2)>log10(pCrit1))==0
                        ii_cross=ii;
                        not_found=0;
                    end
                end
                ii=ii+1;
            end
            
            if ((anal_t_pac(ii_cross)+anal_t_pac(ii_cross-1))/2)-anal_t_pac(ii_t0)>dt
                not_found=1;
            end
            
            if not_found==0
                found_trough_decision_times_pre_control(ii_trough)=1;
                trough_decision_times_pre_control(ii_trough)=((anal_t_pac(ii_cross)+anal_t_pac(ii_cross-1))/2)-t_odor_on;
            end
        else
            
            %Find decision time
            current_ii=find(anal_t_pac>=t_odor_on,1,'first');
            ii_end=current_ii+dii_end;
            
            not_found=1;
            while not_found==1
                delta_ii=find(this_p_val_trough(current_ii:end)<=pCrit,1,'first');
                if current_ii+delta_ii<=ii_end
                    if sum(this_p_val_trough(current_ii+delta_ii:current_ii+delta_ii+dii_crit)>pCrit)==0
                        not_found=0;
                        found_trough_decision_times_pre(ii_trough)=1;
                        trough_decision_times_pre(ii_trough)=anal_t_pac(current_ii+delta_ii)-t_odor_on;
                    else
                        current_ii=current_ii+1;
                        if current_ii>=ii_end
                            not_found=1;
                        end
                    end
                else
                    not_found=0;
                end
            end
            
            %Find control decision time
            current_ii=find(anal_t_pac>=t_odor_on-dt,1,'first');
            ii_end=current_ii+dii_end;
            
            not_found=1;
            while not_found==1
                delta_ii=find(this_p_val_trough(current_ii:end)<=pCrit,1,'first');
                if current_ii+delta_ii<=ii_end
                    if sum(this_p_val_trough(current_ii+delta_ii:current_ii+delta_ii+dii_crit)>pCrit)==0
                        not_found=0;
                        found_trough_decision_times_pre_control(ii_trough)=1;
                        trough_decision_times_pre_control(ii_trough)=anal_t_pac(current_ii+delta_ii)-t_odor_on;
                    else
                        current_ii=current_ii+1;
                        if current_ii>=ii_end
                            not_found=1;
                        end
                    end
                else
                    not_found=0;
                end
            end
        end
    end
    
    fprintf(1, ['For ' bandwidth_names{pacii} ' trough found %d control and %d decision times out of %d mops\n\n'],sum(found_trough_decision_times_pre_control),sum(found_trough_decision_times_pre),length(found_trough_decision_times_pre))
    
    %licks
    these_data=lick_decision_times_hipp(found_lick_decision_times_hipp==1);
    glm_dect.data(glm_dect_ii+1:glm_dect_ii+length(these_data))=these_data;
    glm_dect.lick_peak_trough(glm_dect_ii+1:glm_dect_ii+length(these_data))=zeros(1,length(these_data));
    glm_dect.hipp_vs_pre(glm_dect_ii+1:glm_dect_ii+length(these_data))=ones(1,length(these_data));
    glm_dect.pacii(glm_dect_ii+1:glm_dect_ii+length(these_data))=pacii*ones(1,length(these_data));
    glm_dect_ii=glm_dect_ii+length(these_data);
    
    
    %peak
    these_data=peak_decision_times_pre(found_peak_decision_times_pre==1);
    glm_dect.data(glm_dect_ii+1:glm_dect_ii+length(these_data))=these_data;
    glm_dect.lick_peak_trough(glm_dect_ii+1:glm_dect_ii+length(these_data))=ones(1,length(these_data));
    glm_dect.hipp_vs_pre(glm_dect_ii+1:glm_dect_ii+length(these_data))=ones(1,length(these_data));
    glm_dect.pacii(glm_dect_ii+1:glm_dect_ii+length(these_data))=pacii*ones(1,length(these_data));
    glm_dect_ii=glm_dect_ii+length(these_data);
    
    id_dect_ii=id_dect_ii+1;
    input_dect_data(id_dect_ii).data=these_data;
    input_dect_data(id_dect_ii).description=['Peak prefrontal ' bandwidth_names{pacii}];
    
       %trough
    these_data=trough_decision_times_pre(found_trough_decision_times_pre==1);
    glm_dect.data(glm_dect_ii+1:glm_dect_ii+length(these_data))=these_data;
    glm_dect.lick_peak_trough(glm_dect_ii+1:glm_dect_ii+length(these_data))=2*ones(1,length(these_data));
    glm_dect.hipp_vs_pre(glm_dect_ii+1:glm_dect_ii+length(these_data))=ones(1,length(these_data));
    glm_dect.pacii(glm_dect_ii+1:glm_dect_ii+length(these_data))=pacii*ones(1,length(these_data));
    glm_dect_ii=glm_dect_ii+length(these_data);
    
    id_dect_ii=id_dect_ii+1;
    input_dect_data(id_dect_ii).data=these_data;
    input_dect_data(id_dect_ii).description=['Trough prefrontal ' bandwidth_names{pacii}];
    
    %Now let's plot decision times
    figureNo = figureNo + 1;
    try
        close(figureNo)
    catch
    end
    hFig=figure(figureNo);
    
    ax=gca;ax.LineWidth=3;
    set(hFig, 'units','normalized','position',[.3 .3 .3 .3])
    
    hold on
    
    all_decision_dt_peak=peak_decision_times_pre(found_peak_decision_times_pre==1);
    all_decision_dt_trough=trough_decision_times_pre(found_trough_decision_times_pre==1);
    all_decision_dt_licks=lick_decision_times_hipp(found_lick_decision_times_hipp==1);

    edges=[0:0.033:0.5];
    rand_offset=0.8;
    
    %Peak PRP
    bar_offset=1;
    bar(bar_offset,mean(all_decision_dt_peak),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
    
   
    
    %Trough PRP
    bar_offset=bar_offset+1;
    bar(bar_offset,mean(all_decision_dt_trough),'LineWidth', 3,'EdgeColor','none','FaceColor',[213/255 94/255 0/255])
    

    
    %Licks
    bar_offset=bar_offset+1;
    bar(bar_offset,mean(all_decision_dt_licks),'FaceColor','m','LineWidth', 3,'EdgeColor','none')
    

    %Do the mouse mean calculations and plot the mean points and lines
       
    %licks
    these_data_licks=lick_decision_times_hipp(found_lick_decision_times_hipp==1);
    these_mice_licks=these_mouse_nos_licks(found_lick_decision_times_hipp==1);
    
    %peak
    these_data_peak=peak_decision_times_hipp(found_peak_decision_times_hipp==1);
    these_mice_peak=these_mouse_nos_peak(found_peak_decision_times_hipp==1);
    
    %trough
    these_data_trough=trough_decision_times_hipp(found_trough_decision_times_hipp==1);
    these_mice_trough=these_mouse_nos_trough(found_trough_decision_times_hipp==1);
    
    unique_mouse_nos=unique([these_mice_peak these_mice_trough these_mice_licks]);
    
    for msNo=unique_mouse_nos
        
        %peak
        this_mouse_mean_peak=mean(these_data_peak(these_mice_peak==msNo));
        
        glm_dect_mm.data(glm_ii_mm+1)=this_mouse_mean_peak;
        glm_dect_mm.lick_peak_trough(glm_ii_mm+1)=0;
        glm_dect_mm.hipp_vs_pre(glm_ii_mm+1)=1;
        glm_dect_mm.pacii(glm_ii_mm+1)=pacii;
        glm_ii_mm=glm_ii_mm+1;
        
        %trough
        this_mouse_mean_trough=mean(these_data_trough(these_mice_trough==msNo));
        
        glm_dect_mm.data(glm_ii_mm+1)=this_mouse_mean_trough;
        glm_dect_mm.lick_trough_trough(glm_ii_mm+1)=1;
        glm_dect_mm.hipp_vs_pre(glm_ii_mm+1)=1;
        glm_dect_mm.pacii(glm_ii_mm+1)=pacii;
        glm_ii_mm=glm_ii_mm+1;
        
         %licks
        this_mouse_mean_licks=mean(these_data_licks(these_mice_licks==msNo));
        
        glm_dect_mm.data(glm_ii_mm+1)=this_mouse_mean_licks;
        glm_dect_mm.lick_peak_trough(glm_ii_mm+1)=2;
        glm_dect_mm.hipp_vs_pre(glm_ii_mm+1)=1;
        glm_dect_mm.pacii(glm_ii_mm+1)=pacii;
        glm_ii_mm=glm_ii_mm+1;
        
        
        plot([bar_offset-2 bar_offset-1 bar_offset],[this_mouse_mean_peak this_mouse_mean_trough this_mouse_mean_licks],'-ok','MarkerSize',6,'LineWidth',1,'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7],'Color',[0.7 0.7 0.7])
   
    end
    
    
    %Violin plot
    [mean_out, CIout]=drgViolinPoint(all_decision_dt_peak...
        ,edges,bar_offset-2,rand_offset,'k','k',3);
    
    %Violin plot
    [mean_out, CIout]=drgViolinPoint(all_decision_dt_trough...
        ,edges,bar_offset-1,rand_offset,'k','k',3);
    
    %Violin plot
    [mean_out, CIout]=drgViolinPoint(all_decision_dt_licks...
        ,edges,bar_offset,rand_offset,'k','k',3);
    
    
    title(['Decision time prefrontal ' bandwidth_names{pacii}])
    
    ylabel('dt')
    
    xticks([1 2 3])
    xticklabels({'Peak', 'Trough', 'Licks'})
end



%Perform the glm for decision time
fprintf(1, ['glm for decision time per mouse per odor pair\n'])
fprintf(fileID, ['glm for decision timeper mouse per odor pair\n']);

tbl = table(glm_dect.data',glm_dect.lick_peak_trough',glm_dect.hipp_vs_pre',glm_dect.pacii',...
    'VariableNames',{'detection_time','lick_peak_trough','hipp_vs_pre','bandwidth'});
mdl = fitglm(tbl,'detection_time~lick_peak_trough+hipp_vs_pre+bandwidth+lick_peak_trough*hipp_vs_pre*bandwidth'...
    ,'CategoricalVars',[2,3,4])


txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);

%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for detection time\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values for detection time\n']);

[output_data] = drgMutiRanksumorTtest(input_dect_data, fileID);


%Perform the glm for decision time per mouse
fprintf(1, ['glm for decision time per mouse\n'])
fprintf(fileID, ['glm for decision timeper mouse\n']);
 
tbl = table(glm_dect_mm.data',glm_dect_mm.lick_peak_trough',glm_dect_mm.hipp_vs_pre',glm_dect_mm.pacii',...
    'VariableNames',{'detection_time','lick_peak_trough','hipp_vs_pre','bandwidth'});
mdl = fitglm(tbl,'detection_time~lick_peak_trough+hipp_vs_pre+bandwidth+lick_peak_trough*hipp_vs_pre*bandwidth'...
    ,'CategoricalVars',[2,3,4])


txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


fclose(fileID);

pffft=1;
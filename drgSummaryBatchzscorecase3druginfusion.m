function drgSummaryBatchzscorecase3druginfusion
%Analyzes the z scored PRP analysis performed by
%drgAnalysisBatchLFPCaMKIIv2 case 3
%Takes as a in input the per odorant pair .mat file output for hippocampus or prefrontal

warning('off')

close all hidden
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
grNo=1;

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
PathName='/Users/restrepd/Documents/Projects/CaMKII_analysis/Drug infusion/';

%Text file for statistical output
fileID = fopen([PathName 'drgSummaryBatchzscorecase3DIKN.txt'],'w');

%IMPORTANT: Please note that the numbering of these files is the same as the numbering
%of odorant pairs in the mouse numbering worksheet

%Hippocampus
hippFileName{1}='spm_LFP_druginfusionaceto_one_window_31522_case3__hippoinfLFP80.mat';
hippFileName{2}='spm_LFP_druginfusionethylben_one_window_31522_case3__hippoinfLFP80.mat'; 
% hippFileName{3}='spm_LFP_ethyl_one_win12182021_case3__hippocampusLFP80.mat';
% hippFileName{4}='spm_LFP_PAEA_one_window122121_case3__hippocampusLFP80.mat';
% hippFileName{5}='spm_LFP_pz1ethyl_one_window_122321_case3__hippocampusLFP80.mat';
% hippFileName{6}='spm_LFP_pz1PAEA_one_win12302021_case3__hippocampusLFP80.mat';  
% hippFileName{7}='spm_LFP_pzz1EAPA_one_window01012022_case3__hippocampusLFP80.mat'; 
% hippFileName{8}='spm_LFP_pzz1PAEA_one_window12302021_case3__hippocampusLFP80.mat';
% 
% %Prefrontal
% preFileName{1}='spm_LFP_APEB_one_window_12102021_case3_prefrontalLFP80.mat';
% preFileName{2}='spm_LFP_EBAP_one_window_12192021_case3_prefrontalLFP80.mat'; 
% preFileName{3}='spm_LFP_ethyl_one_win12182021_case3_prefrontalLFP80.mat';
% preFileName{4}='spm_LFP_PAEA_one_window122121_case3_prefrontalLFP80.mat';
% preFileName{5}='spm_LFP_pz1ethyl_one_window_122321_case3_prefrontalLFP80.mat';
% preFileName{6}='spm_LFP_pz1PAEA_one_win12302021_case3_prefrontalLFP80.mat'; 
% preFileName{7}='spm_LFP_pzz1EAPA_one_window01012022_case3_prefrontalLFP80.mat'; 
% preFileName{8}='spm_LFP_pzz1PAEA_one_window12302021_case3_prefrontalLFP80.mat';

 
 
%Load the table of mouse numbers
%Note: This may need to be revised for PRP
mouse_no_table='/Users/restrepd/Documents/Projects/CaMKII_analysis/Reply_to_reviewers/camkii_mice_per_odor_pair_for_DIKN.xlsx';
T_mouse_no = readtable(mouse_no_table);

%Load data hippocampus
all_hippo=[];

for ii=1:length(hippFileName)
    load([PathName hippFileName{ii}])
    all_hippo(ii).handles_out=handles_out;
end

% %Load data prefrontal
% all_pre=[];
% 
% for ii=1:length(preFileName)
%     load([PathName preFileName{ii}])
%     all_pre(ii).handles_out=handles_out;
% end

glm_dect=[];
glm_dect_ii=0;

glm_dect_mm=[];
glm_ii_mm=0;

id_dect_ii=0;
input_dect_data=[];

glm_from=0.5;
glm_to=2.5;



%Now plot the pre p value timecourses for proficient mice for hippocampus
per_ii=1;
anal_t_pac=all_hippo(1).handles_out.anal_t_pac;
for pacii=[1 2]    %for amplitude bandwidths (beta, high gamma)
    
    glm_dect=[];
    glm_dect_ii=0;
    
    glm_dect_mm=[];
    glm_ii_mm=0;
    
    id_dect_ii=0;
    input_dect_data=[];
    
    glm_dect_licks=[];
    glm_dect_licks_ii=0;
    
    id_dect_licks_ii=0;
    input_dect_licks_data=[];
    
    glm_dect_all=[];
    glm_dect_all_ii=0;
    
    id_dect_all_ii=0;
    input_dect_all_data=[];
    
    figureNo = figureNo + 1;
    try
        close(figureNo)
    catch
    end
    hFig=figure(figureNo);
    
    figureNo_peak_trough=figureNo;
    
    ax=gca;ax.LineWidth=3;
    set(hFig, 'units','normalized','position',[.1 .5 .7*0.66 .4])
    
    hold on
    
    
    if pacii==1
        all_lick_dect=[];
        all_lick_found_dect=[];
    end
    all_peak_dect=[];
    all_peak_found_dect=[];
    all_trough_dect=[];
    all_trough_found_dect=[];
    all_dect_groups=[];
    
    if pacii==1
        %This is the lick figure
        figureNo = figureNo + 1;
        try
            close(figureNo)
        catch
        end
        hFig=figure(figureNo);
        
        figureNo_lick=figureNo;
        
        ax=gca;ax.LineWidth=3;
        set(hFig, 'units','normalized','position',[.1 .5 .7*0.33 .4])
        
        hold on
    end
    
    for grNo=1:2
        
        
        %Get these p vals for post-lick
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
                if ~isempty(all_hippo(ii).handles_out.RPtimecourse(these_jjs(jj)).p_vals_peak_above)
                    this_pval=all_hippo(ii).handles_out.RPtimecourse(these_jjs(jj)).p_vals_peak_above;
                    log10_this_pval=zeros(1,length(this_pval));
                    log10_this_pval(1,:)=log10(this_pval);
                    log10_this_pval(isinf(log10_this_pval))=min(log10_this_pval(~isinf(log10_this_pval)));
                    ii_peak=ii_peak+1;
                    all_p_vals_peak(ii_peak,:)=log10_this_pval;
                    this_nn=find(all_hippo(ii).handles_out.RPtimecourse(these_jjs(jj)).mouseNo==these_mouse_no_per_op);
                    these_mouse_nos_peak(ii_peak)=these_mouse_no(this_nn);
                end
                
                %trough
                if ~isempty(all_hippo(ii).handles_out.RPtimecourse(these_jjs(jj)).p_vals_trough_above)
                    this_pval=all_hippo(ii).handles_out.RPtimecourse(these_jjs(jj)).p_vals_trough_above;
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
        peak_decision_times_post_hipp=zeros(1,size(all_p_vals_peak,1));
        found_peak_decision_times_post_hipp=zeros(1,size(all_p_vals_peak,1));
        peak_decision_times_post_hipp_control=zeros(1,size(all_p_vals_peak,1));
        found_peak_decision_times_post_hipp_control=zeros(1,size(all_p_vals_peak,1));
        
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
                    found_peak_decision_times_post_hipp(ii_peak)=1;
                    peak_decision_times_post_hipp(ii_peak)=((anal_t_pac(ii_cross)+anal_t_pac(ii_cross-1))/2)-t_odor_on;
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
                    found_peak_decision_times_post_hipp_control(ii_peak)=1;
                    peak_decision_times_post_hipp_control(ii_peak)=((anal_t_pac(ii_cross)+anal_t_pac(ii_cross-1))/2)-t_odor_on;
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
                            found_peak_decision_times_post_hipp(ii_peak)=1;
                            peak_decision_times_post_hipp(ii_peak)=anal_t_pac(current_ii+delta_ii)-t_odor_on;
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
                            found_peak_decision_times_post_hipp_control(ii_peak)=1;
                            peak_decision_times_post_hipp_control(ii_peak)=anal_t_pac(current_ii+delta_ii)-t_odor_on;
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
        
        fprintf(1, ['For ' bandwidth_names{pacii} ' peak found %d control and %d decision times out of %d mops\n'],sum(found_peak_decision_times_post_hipp_control),sum(found_peak_decision_times_post_hipp),length(found_peak_decision_times_post_hipp))
        
        %trough
        trough_decision_times_post_hipp=zeros(1,size(all_p_vals_trough,1));
        found_trough_decision_times_post_hipp=zeros(1,size(all_p_vals_trough,1));
        trough_decision_times_post_hipp_control=zeros(1,size(all_p_vals_trough,1));
        found_trough_decision_times_post_hipp_control=zeros(1,size(all_p_vals_trough,1));
        
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
                    found_trough_decision_times_post_hipp(ii_trough)=1;
                    trough_decision_times_post_hipp(ii_trough)=((anal_t_pac(ii_cross)+anal_t_pac(ii_cross-1))/2)-t_odor_on;
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
                    found_trough_decision_times_post_hipp_control(ii_trough)=1;
                    trough_decision_times_post_hipp_control(ii_trough)=((anal_t_pac(ii_cross)+anal_t_pac(ii_cross-1))/2)-t_odor_on;
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
                            found_trough_decision_times_post_hipp(ii_trough)=1;
                            trough_decision_times_post_hipp(ii_trough)=anal_t_pac(current_ii+delta_ii)-t_odor_on;
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
                            found_trough_decision_times_post_hipp_control(ii_trough)=1;
                            trough_decision_times_post_hipp_control(ii_trough)=anal_t_pac(current_ii+delta_ii)-t_odor_on;
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
        
        fprintf(1, ['For ' bandwidth_names{pacii} ' trough found %d control and %d decision times out of %d mops\n\n'],sum(found_trough_decision_times_post_hipp_control),sum(found_trough_decision_times_post_hipp),length(found_trough_decision_times_post_hipp))
  
        
        %Save data for glm for decision times
        
        %licks
        these_data=lick_decision_times_hipp(found_lick_decision_times_hipp==1);
        these_mice=these_mouse_nos_licks(found_lick_decision_times_hipp==1);
        glm_dect_licks.data(glm_dect_licks_ii+1:glm_dect_licks_ii+length(these_data))=these_data;
        glm_dect_licks.mouse_no(glm_dect_licks_ii+1:glm_dect_licks_ii+length(these_data))=these_mice;
        glm_dect_licks.lick_peak_trough(glm_dect_licks_ii+1:glm_dect_licks_ii+length(these_data))=zeros(1,length(these_data));
        glm_dect_licks.pre_vs_post(glm_dect_licks_ii+1:glm_dect_licks_ii+length(these_data))=ones(1,length(these_data));
        glm_dect_licks.group(glm_dect_licks_ii+1:glm_dect_licks_ii+length(these_data))=grNo*ones(1,length(these_data));
        glm_dect_licks_ii=glm_dect_licks_ii+length(these_data);
        
        id_dect_licks_ii=id_dect_licks_ii+1;
        input_dect_licks_data(id_dect_licks_ii).data=these_data;
        input_dect_licks_data(id_dect_licks_ii).description=['Licks ' group_legend{grNo}];
        
         %licks
        these_data=lick_decision_times_hipp(found_lick_decision_times_hipp==1);
        these_mice=these_mouse_nos_licks(found_lick_decision_times_hipp==1);
        glm_dect.data(glm_dect_ii+1:glm_dect_ii+length(these_data))=these_data;
        glm_dect.mouse_no(glm_dect_ii+1:glm_dect_ii+length(these_data))=these_mice;
        glm_dect.lick_peak_trough(glm_dect_ii+1:glm_dect_ii+length(these_data))=zeros(1,length(these_data));
        glm_dect.pre_vs_post(glm_dect_ii+1:glm_dect_ii+length(these_data))=ones(1,length(these_data));
        glm_dect.group(glm_dect_ii+1:glm_dect_ii+length(these_data))=grNo*ones(1,length(these_data));
        glm_dect_ii=glm_dect_ii+length(these_data);
        
        id_dect_ii=id_dect_ii+1;
        input_dect_data(id_dect_ii).data=these_data;
        input_dect_data(id_dect_ii).description=['Licks ' group_legend{grNo}];
        
        %peak
        these_data=peak_decision_times_post_hipp(found_peak_decision_times_post_hipp==1);
        these_mice=these_mouse_nos_peak(found_peak_decision_times_post_hipp==1);
        glm_dect.data(glm_dect_ii+1:glm_dect_ii+length(these_data))=these_data;
        glm_dect.mouse_no(glm_dect_ii+1:glm_dect_ii+length(these_data))=these_mice;
        glm_dect.lick_peak_trough(glm_dect_ii+1:glm_dect_ii+length(these_data))=ones(1,length(these_data));
        glm_dect.pre_vs_post(glm_dect_ii+1:glm_dect_ii+length(these_data))=ones(1,length(these_data));
        glm_dect.group(glm_dect_ii+1:glm_dect_ii+length(these_data))=grNo*ones(1,length(these_data));
        glm_dect_ii=glm_dect_ii+length(these_data);
        
        id_dect_ii=id_dect_ii+1;
        input_dect_data(id_dect_ii).data=these_data;
        input_dect_data(id_dect_ii).description=['Peak post-lick ' group_legend{grNo}];
        
        %trough
        these_data=trough_decision_times_post_hipp(found_trough_decision_times_post_hipp==1);
        these_mice=these_mouse_nos_trough(found_trough_decision_times_post_hipp==1);
        glm_dect.data(glm_dect_ii+1:glm_dect_ii+length(these_data))=these_data;
        glm_dect.mouse_no(glm_dect_ii+1:glm_dect_ii+length(these_data))=these_mice;
        glm_dect.lick_peak_trough(glm_dect_ii+1:glm_dect_ii+length(these_data))=2*ones(1,length(these_data));
        glm_dect.pre_vs_post(glm_dect_ii+1:glm_dect_ii+length(these_data))=ones(1,length(these_data));
        glm_dect.group(glm_dect_ii+1:glm_dect_ii+length(these_data))=grNo*ones(1,length(these_data));
        glm_dect_ii=glm_dect_ii+length(these_data);
        
        id_dect_ii=id_dect_ii+1;
        input_dect_data(id_dect_ii).data=these_data;
        input_dect_data(id_dect_ii).description=['Trough post-lick ' group_legend{grNo}];
        
        %Now do glm_dect_all
         
        
        %peak
        these_data=peak_decision_times_post_hipp(found_peak_decision_times_post_hipp==1);
        these_mice=these_mouse_nos_peak(found_peak_decision_times_post_hipp==1);
        glm_dect_all.data(glm_dect_all_ii+1:glm_dect_all_ii+length(these_data))=these_data;
        glm_dect_all.mouse_no(glm_dect_all_ii+1:glm_dect_all_ii+length(these_data))=these_mice;
        glm_dect_all.lick_peak_trough(glm_dect_all_ii+1:glm_dect_all_ii+length(these_data))=ones(1,length(these_data));
        glm_dect_all.pre_vs_post(glm_dect_all_ii+1:glm_dect_all_ii+length(these_data))=ones(1,length(these_data));
        glm_dect_all.group(glm_dect_all_ii+1:glm_dect_all_ii+length(these_data))=grNo*ones(1,length(these_data));
        glm_dect_all.hipp_vs_pre(glm_dect_all_ii+1:glm_dect_all_ii+length(these_data))=zeros(1,length(these_data));
        glm_dect_all_ii=glm_dect_all_ii+length(these_data);
        
        id_dect_all_ii=id_dect_all_ii+1;
        input_dect_all_data(id_dect_all_ii).data=these_data;
        input_dect_all_data(id_dect_all_ii).description=['Peak post-lick ' group_legend{grNo} ' hippocampus'];
        
        %trough
        these_data=trough_decision_times_post_hipp(found_trough_decision_times_post_hipp==1);
        these_mice=these_mouse_nos_trough(found_trough_decision_times_post_hipp==1);
        glm_dect_all.data(glm_dect_all_ii+1:glm_dect_all_ii+length(these_data))=these_data;
        glm_dect_all.mouse_no(glm_dect_all_ii+1:glm_dect_all_ii+length(these_data))=these_mice;
        glm_dect_all.lick_peak_trough(glm_dect_all_ii+1:glm_dect_all_ii+length(these_data))=2*ones(1,length(these_data));
        glm_dect_all.pre_vs_post(glm_dect_all_ii+1:glm_dect_all_ii+length(these_data))=ones(1,length(these_data));
        glm_dect_all.group(glm_dect_all_ii+1:glm_dect_all_ii+length(these_data))=grNo*ones(1,length(these_data));
        glm_dect_all.hipp_vs_pre(glm_dect_all_ii+1:glm_dect_all_ii+length(these_data))=zeros(1,length(these_data));
        glm_dect_all_ii=glm_dect_all_ii+length(these_data);
        
        id_dect_all_ii=id_dect_all_ii+1;
        input_dect_all_data(id_dect_all_ii).data=these_data;
        input_dect_all_data(id_dect_all_ii).description=['Trough post-lick ' group_legend{grNo} ' hippocampus'];
        
        
        %Now let's plot decision times
        all_decision_dt_peak=peak_decision_times_post_hipp(found_peak_decision_times_post_hipp==1);
        all_decision_dt_trough=trough_decision_times_post_hipp(found_trough_decision_times_post_hipp==1);
        all_decision_dt_licks=lick_decision_times_hipp(found_lick_decision_times_hipp==1);
        
        edges=[0:0.033:0.5];
        rand_offset=0.8;
        
        %Peak PRP
        figure(figureNo_peak_trough)
        hold on
        bar_offset=grNo;
        switch grNo
            case 1
                bar(bar_offset,mean(all_decision_dt_peak),'g','LineWidth', 3,'EdgeColor','none')
            case 2
                bar(bar_offset,mean(all_decision_dt_peak),'b','LineWidth', 3,'EdgeColor','none')
            case 3
                bar(bar_offset,mean(all_decision_dt_peak),'y','LineWidth', 3,'EdgeColor','none')
        end

        %Trough PRP
        bar_offset=grNo+4;
        switch grNo
            case 1
                bar(bar_offset,mean(all_decision_dt_trough),'g','LineWidth', 3,'EdgeColor','none')
            case 2
                bar(bar_offset,mean(all_decision_dt_trough),'b','LineWidth', 3,'EdgeColor','none')
            case 3
                bar(bar_offset,mean(all_decision_dt_trough),'y','LineWidth', 3,'EdgeColor','none')
        end
        
        
        if pacii==1
            %Licks
            figure(figureNo_lick)
            hold on
            bar_offset=grNo+8;
            switch grNo
                case 1
                    bar(bar_offset,mean(all_decision_dt_licks),'g','LineWidth', 3,'EdgeColor','none')
                case 2
                    bar(bar_offset,mean(all_decision_dt_licks),'b','LineWidth', 3,'EdgeColor','none')
                case 3
                    bar(bar_offset,mean(all_decision_dt_licks),'y','LineWidth', 3,'EdgeColor','none')
            end
        end
        
        %Do the mouse mean calculations and plot the mean points and lines
        
        %licks
        these_data_licks=lick_decision_times_hipp(found_lick_decision_times_hipp==1);
        these_mice_licks=these_mouse_nos_licks(found_lick_decision_times_hipp==1);
        
        %peak
        these_data_peak=peak_decision_times_post_hipp(found_peak_decision_times_post_hipp==1);
        these_mice_peak=these_mouse_nos_peak(found_peak_decision_times_post_hipp==1);
        
        %trough
        these_data_trough=trough_decision_times_post_hipp(found_trough_decision_times_post_hipp==1);
        these_mice_trough=these_mouse_nos_trough(found_trough_decision_times_post_hipp==1);
        
        unique_mouse_nos=unique([these_mice_peak these_mice_trough these_mice_licks]);
        
        %Save data for correlation plots
        all_peak_dect=[all_peak_dect peak_decision_times_post_hipp];
        all_peak_found_dect=[all_peak_found_dect found_peak_decision_times_post_hipp];
        
        all_trough_dect=[all_trough_dect trough_decision_times_post_hipp];
        all_trough_found_dect=[all_trough_found_dect found_trough_decision_times_post_hipp];
        
        all_dect_groups=[all_dect_groups grNo*ones(1,length(found_peak_decision_times_post_hipp))];
        
        if pacii==1
            all_lick_dect=[all_lick_dect lick_decision_times_hipp];
            all_lick_found_dect=[all_lick_found_dect found_lick_decision_times_hipp];
        end
        
        fprintf(1, ['Number of mice for ' group_legend{grNo} ' = ' num2str(length(unique_mouse_nos)) '\n\n']);
        fprintf(fileID, ['Number of mice for ' group_legend{grNo} ' = ' num2str(length(unique_mouse_nos)) '\n\n']);
        
        for msNo=unique_mouse_nos
            
            %peak
            this_mouse_mean_peak=mean(these_data_peak(these_mice_peak==msNo));
            
            glm_dect_mm.data(glm_ii_mm+1)=this_mouse_mean_peak;
            glm_dect_mm.lick_peak_trough(glm_ii_mm+1)=0;
            glm_dect_mm.pre_vs_post(glm_ii_mm+1)=1;
            glm_ii_mm=glm_ii_mm+1;
            
            %trough
            this_mouse_mean_trough=mean(these_data_trough(these_mice_trough==msNo));
            
            glm_dect_mm.data(glm_ii_mm+1)=this_mouse_mean_trough;
            glm_dect_mm.lick_peak_trough(glm_ii_mm+1)=1;
            glm_dect_mm.pre_vs_post(glm_ii_mm+1)=1;
            glm_ii_mm=glm_ii_mm+1;
            
            %licks
            this_mouse_mean_licks=mean(these_data_licks(these_mice_licks==msNo));
            
            glm_dect_mm.data(glm_ii_mm+1)=this_mouse_mean_licks;
            glm_dect_mm.lick_peak_trough(glm_ii_mm+1)=2;
            glm_dect_mm.pre_vs_post(glm_ii_mm+1)=1;
            glm_ii_mm=glm_ii_mm+1;
            
            figure(figureNo_peak_trough)
            hold on
            plot([grNo grNo+4],[this_mouse_mean_peak this_mouse_mean_trough],'o','MarkerSize',6,'MarkerFaceColor',[0.6 0.6 0.6],'MarkerEdgeColor',[0.6 0.6 0.6])
            
            if pacii==1
                figure(figureNo_lick)
                hold on
                plot([grNo+8],[this_mouse_mean_licks],'o','MarkerSize',6,'MarkerFaceColor',[0.6 0.6 0.6],'MarkerEdgeColor',[0.6 0.6 0.6])
            end
            
        end
        
        figure(figureNo_peak_trough)
            hold on
            
        %Violin plot
        [mean_out, CIout]=drgViolinPoint(all_decision_dt_peak...
            ,edges,grNo,rand_offset,'k','k',3);
        
        
         
        %Violin plot
        try
        [mean_out, CIout]=drgViolinPoint(all_decision_dt_trough...
            ,edges,grNo+4,rand_offset,'k','k',3);
        catch
        end
        
        
        if pacii==1
            figure(figureNo_lick)
            hold on
            
            %Violin plot
            [mean_out, CIout]=drgViolinPoint(all_decision_dt_licks...
                ,edges,grNo+8,rand_offset,'k','k',3);
            
        end
        
       

    end
    
   
    figure(figureNo_peak_trough)
    hold on
    
    title(['Decision time post-lick ' bandwidth_names{pacii} ' hippocampus'])
    ylim([0 0.8])
    ylabel('dt')
    
    xticks([2 6])
    xticklabels({'Peak', 'Trough'})
    
    if pacii==1
        figure(figureNo_lick)
        hold on
        
        title(['Decision time for licks ' bandwidth_names{pacii} ' hippocampus'])
        ylim([0 0.8])
        ylabel('dt')
        
        
        
        %Perform the glm for decision time for licks
        fprintf(1, ['glm for lick decision time per mouse per odor pair \n']);
        fprintf(fileID, ['glm for lick decision time per mouse per odor pair \n']);
        
        tbl = table(glm_dect_licks.data',glm_dect_licks.group',...
            'VariableNames',{'detection_time','genotype'});
        mdl = fitglm(tbl,'detection_time~genotype'...
            ,'CategoricalVars',[2])
        
        
        txt = evalc('mdl');
        txt=regexp(txt,'<strong>','split');
        txt=cell2mat(txt);
        txt=regexp(txt,'</strong>','split');
        txt=cell2mat(txt);
        
        fprintf(fileID,'%s\n', txt);
        
        %Do the ranksum/t-test
        fprintf(1, ['\n\nRanksum or t-test p values for lick decision time per mouse per odor pair \n']);
        fprintf(fileID, ['\n\nRanksum or t-test p values for lick decision time per mouse per odor pair \n']);
        
        
        [output_data] = drgMutiRanksumorTtest(input_dect_licks_data, fileID);
        
        %Nested ANOVAN
        %https://www.mathworks.com/matlabcentral/answers/491365-within-between-subjects-in-anovan
        nesting=[0 0 ; ... % This line indicates that group factor is not nested in any other factor.
            1 0];    % This line indicates that mouse_no (the last factor) is nested under group, perCorr and event
        % (the 1 in position 1 on the line indicates nesting under the first factor).
        figureNo = figureNo + 1;
        
        [p anovanTbl stats]=anovan(glm_dect_licks.data,{glm_dect_licks.group glm_dect_licks.mouse_no},...
            'model','interaction',...
            'nested',nesting,...
            'varnames',{'genotype','mouse_no'});
        
        fprintf(fileID, ['\n\nNested ANOVAN for lick decision time per mouse per odor pair\n']);
        drgWriteANOVANtbl(anovanTbl,fileID);
    end
    
    %Perform the glm for decision time
    fprintf(1, ['glm for post-lick decision time per mouse per odor pair for hippocampus ' bandwidth_names{pacii} '\n']);
    fprintf(fileID, ['glm for post-lick decision time per mouse per odor pair for hippocampus ' bandwidth_names{pacii} '\n']);
    
    tbl = table(glm_dect.data',glm_dect.lick_peak_trough',glm_dect.group',...
        'VariableNames',{'detection_time','lick_peak_trough','genotype'});
    mdl = fitglm(tbl,'detection_time~genotype+lick_peak_trough+lick_peak_trough*genotype'...
        ,'CategoricalVars',[2,3])
    
    
    txt = evalc('mdl');
    txt=regexp(txt,'<strong>','split');
    txt=cell2mat(txt);
    txt=regexp(txt,'</strong>','split');
    txt=cell2mat(txt);
    
    fprintf(fileID,'%s\n', txt);
    
    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for post-lick decision time per mouse per odor pair for hippocampus '  bandwidth_names{pacii} '\n']);
    fprintf(fileID, ['\n\nRanksum or t-test p values for post-lick decision time per mouse per odor pair for hippocampus '  bandwidth_names{pacii} '\n']);
     
    
    [output_data] = drgMutiRanksumorTtest(input_dect_data, fileID);
    
    %Nested ANOVAN
    %https://www.mathworks.com/matlabcentral/answers/491365-within-between-subjects-in-anovan
    nesting=[0 0 0 ; ... % This line indicates that group factor is not nested in any other factor.
        0 0 0 ; ... % This line indicates that perCorr is not nested in any other factor.
        1 1 0];    % This line indicates that mouse_no (the last factor) is nested under group, perCorr and event
    % (the 1 in position 1 on the line indicates nesting under the first factor).
    figureNo = figureNo + 1;
    
    [p anovanTbl stats]=anovan(glm_dect.data,{glm_dect.lick_peak_trough glm_dect.group glm_dect.mouse_no},...
        'model','interaction',...
        'nested',nesting,...
        'varnames',{'lick_peak_trough','genotype','mouse_no'});
    
    fprintf(fileID, ['\n\nNested ANOVAN for post-lick decision time per mouse per odor pair for '  bandwidth_names{pacii} ' hippocampus\n']);
    drgWriteANOVANtbl(anovanTbl,fileID);
    
    %Plot the correlations for the peak
    figureNo = figureNo + 1;
    try
        close(figureNo)
    catch
    end
    hFig=figure(figureNo);
    
    set(hFig, 'units','normalized','position',[.1 .5 .25 .4])
    
    hold on
    
    ax=gca;ax.LineWidth=3;
    
    for grNo=1:2
        these_included=all_peak_found_dect&all_lick_found_dect&(all_dect_groups==grNo);
        switch grNo
            case 1
                plot(all_peak_dect(these_included),all_lick_dect(these_included),'o','MarkerEdgeColor','g','MarkerFaceColor','g')
            case 2
                plot(all_peak_dect(these_included),all_lick_dect(these_included),'o','MarkerEdgeColor','b','MarkerFaceColor','b')
            case 3
                plot(all_peak_dect(these_included),all_lick_dect(these_included),'o','MarkerEdgeColor','y','MarkerFaceColor','y')
        end
        
    end
    
    these_included=all_peak_found_dect&all_lick_found_dect;
    y=all_lick_dect(these_included)';
    x=all_peak_dect(these_included)';
    
    
    [rho,pval]=corr(x,y);
    
    %Fit a line
    c = polyfit(x,y,1);
    % Display evaluated equation y = m*x + b
    disp(['Equation is y = ' num2str(c(1)) '*x + ' num2str(c(2))])
    % Evaluate fit equation using polyval
    y_est = polyval(c,x);
    % Add trend line to plot
    plot(x,y_est,'k-','LineWidth',2)
    
    xlim([0 0.8])
    ylim([0 0.8])
    
    fprintf(1, ['\nrho = %d, p value = %d for hippocampus peak post-lick ztPRP vs lick discrimination time for '  bandwidth_names{pacii} '\n\n'],rho,pval)
    
    ylabel('Discrimintion time for licks')
    xlabel('Discrimination time for peak pre-lick ztPRP')
    title(['Discrimination time for licks vs. peak post-lick ztPRP '  bandwidth_names{pacii} ' hippocampus' ])
    
    %Plot the correlations for the trough
    figureNo = figureNo + 1;
    try
        close(figureNo)
    catch
    end
    hFig=figure(figureNo);
    
    set(hFig, 'units','normalized','position',[.1 .5 .25 .4])
    
    hold on
    
    ax=gca;ax.LineWidth=3;
    
    for grNo=1:2
        these_included=all_trough_found_dect&all_lick_found_dect&(all_dect_groups==grNo);
        switch grNo
            case 1
                plot(all_trough_dect(these_included),all_lick_dect(these_included),'o','MarkerEdgeColor','g','MarkerFaceColor','g')
            case 2
                plot(all_trough_dect(these_included),all_lick_dect(these_included),'o','MarkerEdgeColor','b','MarkerFaceColor','b')
            case 3
                plot(all_trough_dect(these_included),all_lick_dect(these_included),'o','MarkerEdgeColor','y','MarkerFaceColor','y')
        end
        
    end
    
    these_included=all_trough_found_dect&all_lick_found_dect;
    y=all_lick_dect(these_included)';
    x=all_trough_dect(these_included)';
    
    
    [rho,pval]=corr(x,y);
    
    %Fit a line
    c = polyfit(x,y,1);
    % Display evaluated equation y = m*x + b
    disp(['Equation is y = ' num2str(c(1)) '*x + ' num2str(c(2))])
    % Evaluate fit equation using polyval
    y_est = polyval(c,x);
    % Add trend line to plot
    plot(x,y_est,'k-','LineWidth',2)
    
    xlim([0 0.8])
    ylim([0 0.8])
    
    fprintf(1, ['\nrho = %d, p value = %d for hippocampus trough post-lick ztPRP vs lick discrimination time for '  bandwidth_names{pacii} '\n\n'],rho,pval)
    
    ylabel('Discrimintion time for licks')
    xlabel('Discrimination time for trough pre-lick ztPRP')
    title(['Discrimination time for trough post-lick ztPRP vs licks '  bandwidth_names{pacii} ' hippocampus' ])
    
    
    %Now do prefrontal
    glm_dect=[];
    glm_dect_ii=0;
    
    glm_dect_mm=[];
    glm_ii_mm=0;
    
    id_dect_ii=0;
    input_dect_data=[];
    
    figureNo = figureNo + 1;
    try
        close(figureNo)
    catch
    end
    hFig=figure(figureNo);
    
    ax=gca;ax.LineWidth=3;
    set(hFig, 'units','normalized','position',[.1 .5 .7*0.66 .4])
    
    hold on
    
      
    if pacii==1
        all_lick_dect=[];
        all_lick_found_dect=[];
    end
    all_peak_dect=[];
    all_peak_found_dect=[];
    all_trough_dect=[];
    all_trough_found_dect=[];
    all_dect_groups=[];
    
    for grNo=1:2
    
        
        
        %Get these p vals
        all_p_vals_peak_above=[];
        ii_peak=0;
        all_p_vals_trough_above=[];
        ii_trough=0;
        all_p_vals_licks=[];
        ii_licks=0;
         
        for ii=1:length(preFileName)
            these_jjs=[];
            these_mouse_no_per_op=T_mouse_no.mouse_no_per_op(T_mouse_no.odor_pair_no==ii);
            these_mouse_no=T_mouse_no.mouse_no(T_mouse_no.odor_pair_no==ii);
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
                if ~isempty(all_pre(ii).handles_out.RPtimecourse(these_jjs(jj)).p_vals_peak_above)
                    this_pval=all_pre(ii).handles_out.RPtimecourse(these_jjs(jj)).p_vals_peak_above;
                    log10_this_pval=zeros(1,length(this_pval));
                    log10_this_pval(1,:)=log10(this_pval);
                    log10_this_pval(isinf(log10_this_pval))=min(log10_this_pval(~isinf(log10_this_pval)));
                    ii_peak=ii_peak+1;
                    all_p_vals_peak_above(ii_peak,:)=log10_this_pval;
                    this_nn=find(all_pre(ii).handles_out.RPtimecourse(these_jjs(jj)).mouseNo==these_mouse_no_per_op);
                    these_mouse_nos_peak(ii_peak)=these_mouse_no(this_nn);
                end
                
                %trough
                if ~isempty(all_pre(ii).handles_out.RPtimecourse(these_jjs(jj)).p_vals_trough_above)
                    this_pval=all_pre(ii).handles_out.RPtimecourse(these_jjs(jj)).p_vals_trough_above;
                    log10_this_pval=zeros(1,length(this_pval));
                    log10_this_pval(1,:)=log10(this_pval);
                    log10_this_pval(isinf(log10_this_pval))=min(log10_this_pval(~isinf(log10_this_pval)));
                    ii_trough=ii_trough+1;
                    all_p_vals_trough_above(ii_trough,:)=log10_this_pval;
                    this_nn=find(all_pre(ii).handles_out.RPtimecourse(these_jjs(jj)).mouseNo==these_mouse_no_per_op);
                    these_mouse_nos_trough(ii_peak)=these_mouse_no(this_nn);
                end
                
                %licks
                if ~isempty(all_pre(ii).handles_out.RPtimecourse(these_jjs(jj)).p_vals_licks)
                    this_pval=all_pre(ii).handles_out.RPtimecourse(these_jjs(jj)).p_vals_licks;
                    if sum(isnan(this_pval))~=0
                        ii_nan=find(isnan(this_pval));
                        for kk=ii_nan
                            this_pval(kk)=this_pval(kk-1);
                        end
                    end
                    log10_this_pval=zeros(1,length(this_pval));
                    log10_this_pval(1,:)=log10(this_pval);
                    log10_this_pval(isinf(log10_this_pval))=min(log10_this_pval(~isinf(log10_this_pval)));
                    ii_licks=ii_licks+1;
                    all_p_vals_licks(ii_licks,:)=log10_this_pval;
                    this_nn=find(all_pre(ii).handles_out.RPtimecourse(these_jjs(jj)).mouseNo==these_mouse_no_per_op);
                    these_mouse_nos_licks(ii_peak)=these_mouse_no(this_nn);
                end
                
            end
        end
        
        %Find the decision making times
        dii_crit=floor(dt_crit/(anal_t_pac(2)-anal_t_pac(1)));
        dii_end=floor(dt/(anal_t_pac(2)-anal_t_pac(1)));
        
        lick_decision_times_pre=zeros(1,size(all_p_vals_licks,1));
        found_lick_decision_times_pre=zeros(1,size(all_p_vals_licks,1));
        lick_decision_times_pre_control=zeros(1,size(all_p_vals_licks,1));
        found_lick_decision_times_pre_control=zeros(1,size(all_p_vals_licks,1));
        
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
                    found_lick_decision_times_pre(ii_lick)=1;
                    lick_decision_times_pre(ii_lick)=((anal_t_pac(ii_cross)+anal_t_pac(ii_cross-1))/2)-t_odor_on;
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
                    found_lick_decision_times_pre_control(ii_lick)=1;
                    lick_decision_times_pre_control(ii_lick)=((anal_t_pac(ii_cross)+anal_t_pac(ii_cross-1))/2)-t_odor_on;
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
                            found_lick_decision_times_pre(ii_lick)=1;
                            lick_decision_times_pre(ii_lick)=anal_t_pac(current_ii+delta_ii)-t_odor_on;
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
                            found_lick_decision_times_pre_control(ii_lick)=1;
                            lick_decision_times_pre_control(ii_lick)=anal_t_pac(current_ii+delta_ii)-t_odor_on;
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
        
        fprintf(1, ['For ' bandwidth_names{pacii} ' licks found %d control and %d decision times out of %d mops\n'],sum(found_lick_decision_times_pre_control),sum(found_lick_decision_times_pre),length(found_lick_decision_times_pre))
        
        %peak
        peak_decision_times_post_pref=zeros(1,size(all_p_vals_peak_above,1));
        found_peak_decision_times_post_pref=zeros(1,size(all_p_vals_peak_above,1));
        peak_decision_times_post_pref_control=zeros(1,size(all_p_vals_peak_above,1));
        found_peak_decision_times_post_pref_control=zeros(1,size(all_p_vals_peak_above,1));
        
        for ii_peak=1:size(all_p_vals_peak_above,1)
            
            this_p_val_peak=all_p_vals_peak_above(ii_peak,:);
            
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
                    found_peak_decision_times_post_pref(ii_peak)=1;
                    peak_decision_times_post_pref(ii_peak)=((anal_t_pac(ii_cross)+anal_t_pac(ii_cross-1))/2)-t_odor_on;
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
                    found_peak_decision_times_post_pref_control(ii_peak)=1;
                    peak_decision_times_post_pref_control(ii_peak)=((anal_t_pac(ii_cross)+anal_t_pac(ii_cross-1))/2)-t_odor_on;
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
                            found_peak_decision_times_post_pref(ii_peak)=1;
                            peak_decision_times_post_pref(ii_peak)=anal_t_pac(current_ii+delta_ii)-t_odor_on;
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
                            found_peak_decision_times_post_pref_control(ii_peak)=1;
                            peak_decision_times_post_pref_control(ii_peak)=anal_t_pac(current_ii+delta_ii)-t_odor_on;
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
        
        fprintf(1, ['For ' bandwidth_names{pacii} ' peak found %d control and %d decision times out of %d mops\n'],sum(found_peak_decision_times_post_pref_control),sum(found_peak_decision_times_post_pref),length(found_peak_decision_times_post_pref))
        
        %trough
        trough_decision_times_post_pref=zeros(1,size(all_p_vals_trough_above,1));
        found_trough_decision_times_post_pref=zeros(1,size(all_p_vals_trough_above,1));
        trough_decision_times_post_pref_control=zeros(1,size(all_p_vals_trough_above,1));
        found_trough_decision_times_post_pref_control=zeros(1,size(all_p_vals_trough_above,1));
        
        for ii_trough=1:size(all_p_vals_trough_above,1)
            
            this_p_val_trough=all_p_vals_trough_above(ii_trough,:);
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
                    found_trough_decision_times_post_pref(ii_trough)=1;
                    trough_decision_times_post_pref(ii_trough)=((anal_t_pac(ii_cross)+anal_t_pac(ii_cross-1))/2)-t_odor_on;
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
                    found_trough_decision_times_post_pref_control(ii_trough)=1;
                    trough_decision_times_post_pref_control(ii_trough)=((anal_t_pac(ii_cross)+anal_t_pac(ii_cross-1))/2)-t_odor_on;
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
                            found_trough_decision_times_post_pref(ii_trough)=1;
                            trough_decision_times_post_pref(ii_trough)=anal_t_pac(current_ii+delta_ii)-t_odor_on;
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
                            found_trough_decision_times_post_pref_control(ii_trough)=1;
                            trough_decision_times_post_pref_control(ii_trough)=anal_t_pac(current_ii+delta_ii)-t_odor_on;
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
        
        fprintf(1, ['For ' bandwidth_names{pacii} ' trough found %d control and %d decision times out of %d mops\n\n'],sum(found_trough_decision_times_post_pref_control),sum(found_trough_decision_times_post_pref),length(found_trough_decision_times_post_pref))
        
        %licks
        these_data=lick_decision_times_pre(found_lick_decision_times_pre==1);
        these_mice=these_mouse_nos_licks(found_lick_decision_times_pre==1);
        glm_dect.data(glm_dect_ii+1:glm_dect_ii+length(these_data))=these_data;
        glm_dect.mouse_no(glm_dect_ii+1:glm_dect_ii+length(these_data))=these_mice;
        glm_dect.lick_peak_trough(glm_dect_ii+1:glm_dect_ii+length(these_data))=zeros(1,length(these_data));
        glm_dect.pre_vs_post(glm_dect_ii+1:glm_dect_ii+length(these_data))=ones(1,length(these_data));
        glm_dect.group(glm_dect_ii+1:glm_dect_ii+length(these_data))=grNo*ones(1,length(these_data));
        glm_dect_ii=glm_dect_ii+length(these_data);
        
        id_dect_ii=id_dect_ii+1;
        input_dect_data(id_dect_ii).data=these_data;
        input_dect_data(id_dect_ii).description=['Lick post-lick ' group_legend{grNo}];
        
        %peak
        these_data=peak_decision_times_post_pref(found_peak_decision_times_post_pref==1);
        these_mice=these_mouse_nos_peak(found_peak_decision_times_post_pref==1);
        glm_dect.data(glm_dect_ii+1:glm_dect_ii+length(these_data))=these_data;
        glm_dect.mouse_no(glm_dect_ii+1:glm_dect_ii+length(these_data))=these_mice;
        glm_dect.lick_peak_trough(glm_dect_ii+1:glm_dect_ii+length(these_data))=ones(1,length(these_data));
        glm_dect.pre_vs_post(glm_dect_ii+1:glm_dect_ii+length(these_data))=ones(1,length(these_data));
        glm_dect.group(glm_dect_ii+1:glm_dect_ii+length(these_data))=grNo*ones(1,length(these_data));
        glm_dect_ii=glm_dect_ii+length(these_data);
        
        id_dect_ii=id_dect_ii+1;
        input_dect_data(id_dect_ii).data=these_data;
        input_dect_data(id_dect_ii).description=['Peak post-lick ' group_legend{grNo}];
        
        %trough
        these_data=trough_decision_times_post_pref(found_trough_decision_times_post_pref==1);
        these_mice=these_mouse_nos_trough(found_trough_decision_times_post_pref==1);
        glm_dect.data(glm_dect_ii+1:glm_dect_ii+length(these_data))=these_data;
        glm_dect.mouse_no(glm_dect_ii+1:glm_dect_ii+length(these_data))=these_mice;
        glm_dect.lick_peak_trough(glm_dect_ii+1:glm_dect_ii+length(these_data))=2*ones(1,length(these_data));
        glm_dect.pre_vs_post(glm_dect_ii+1:glm_dect_ii+length(these_data))=ones(1,length(these_data));
        glm_dect.group(glm_dect_ii+1:glm_dect_ii+length(these_data))=grNo*ones(1,length(these_data));
        glm_dect_ii=glm_dect_ii+length(these_data);
        
        id_dect_ii=id_dect_ii+1;
        input_dect_data(id_dect_ii).data=these_data;
        input_dect_data(id_dect_ii).description=['Trough post-lick ' group_legend{grNo}];
        
        %Now do glm_dect_all
         
        %peak
        these_data=peak_decision_times_post_pref(found_peak_decision_times_post_pref==1);
        these_mice=these_mouse_nos_peak(found_peak_decision_times_post_pref==1);
        glm_dect_all.data(glm_dect_all_ii+1:glm_dect_all_ii+length(these_data))=these_data;
        glm_dect_all.mouse_no(glm_dect_all_ii+1:glm_dect_all_ii+length(these_data))=these_mice;
        glm_dect_all.lick_peak_trough(glm_dect_all_ii+1:glm_dect_all_ii+length(these_data))=ones(1,length(these_data));
        glm_dect_all.pre_vs_post(glm_dect_all_ii+1:glm_dect_all_ii+length(these_data))=ones(1,length(these_data));
        glm_dect_all.group(glm_dect_all_ii+1:glm_dect_all_ii+length(these_data))=grNo*ones(1,length(these_data));
        glm_dect_all.hipp_vs_pre(glm_dect_all_ii+1:glm_dect_all_ii+length(these_data))=ones(1,length(these_data));
        glm_dect_all_ii=glm_dect_all_ii+length(these_data);
        
        id_dect_all_ii=id_dect_all_ii+1;
        input_dect_all_data(id_dect_all_ii).data=these_data;
        input_dect_all_data(id_dect_all_ii).description=['Peak post-lick ' group_legend{grNo} ' prefrontal'];
        
        %trough
        these_data=trough_decision_times_post_pref(found_trough_decision_times_post_pref==1);
        these_mice=these_mouse_nos_trough(found_trough_decision_times_post_pref==1);
        glm_dect_all.data(glm_dect_all_ii+1:glm_dect_all_ii+length(these_data))=these_data;
        glm_dect_all.mouse_no(glm_dect_all_ii+1:glm_dect_all_ii+length(these_data))=these_mice;
        glm_dect_all.lick_peak_trough(glm_dect_all_ii+1:glm_dect_all_ii+length(these_data))=2*ones(1,length(these_data));
        glm_dect_all.pre_vs_post(glm_dect_all_ii+1:glm_dect_all_ii+length(these_data))=ones(1,length(these_data));
        glm_dect_all.group(glm_dect_all_ii+1:glm_dect_all_ii+length(these_data))=grNo*ones(1,length(these_data));
        glm_dect_all.hipp_vs_pre(glm_dect_all_ii+1:glm_dect_all_ii+length(these_data))=ones(1,length(these_data));
        glm_dect_all_ii=glm_dect_all_ii+length(these_data);
        
        id_dect_all_ii=id_dect_all_ii+1;
        input_dect_all_data(id_dect_all_ii).data=these_data;
        input_dect_all_data(id_dect_all_ii).description=['Trough post-lick ' group_legend{grNo} ' prefrontal'];
        
        %Now let's plot decision times
      
        all_decision_dt_peak=peak_decision_times_post_pref(found_peak_decision_times_post_pref==1);
        all_decision_dt_trough=trough_decision_times_post_pref(found_trough_decision_times_post_pref==1);
        all_decision_dt_licks=lick_decision_times_pre(found_lick_decision_times_pre==1);
        
        edges=[0:0.033:0.5];
        rand_offset=0.8;
        
        
        %Peak PRP
        bar_offset=grNo;
        switch grNo
            case 1
                bar(bar_offset,mean(all_decision_dt_peak),'g','LineWidth', 3,'EdgeColor','none')
            case 2
                bar(bar_offset,mean(all_decision_dt_peak),'b','LineWidth', 3,'EdgeColor','none')
            case 3
                bar(bar_offset,mean(all_decision_dt_peak),'y','LineWidth', 3,'EdgeColor','none')
        end

        %Trough PRP
        bar_offset=grNo+4;
        switch grNo
            case 1
                bar(bar_offset,mean(all_decision_dt_trough),'g','LineWidth', 3,'EdgeColor','none')
            case 2
                bar(bar_offset,mean(all_decision_dt_trough),'b','LineWidth', 3,'EdgeColor','none')
            case 3
                bar(bar_offset,mean(all_decision_dt_trough),'y','LineWidth', 3,'EdgeColor','none')
        end
        
        
        
%         %Licks
%          bar_offset=grNo+8;
%         switch grNo
%             case 1
%                 bar(bar_offset,mean(all_decision_dt_licks),'g','LineWidth', 3,'EdgeColor','none')
%             case 2
%                 bar(bar_offset,mean(all_decision_dt_licks),'b','LineWidth', 3,'EdgeColor','none')
%             case 3
%                 bar(bar_offset,mean(all_decision_dt_licks),'y','LineWidth', 3,'EdgeColor','none')
%         end
        
        %Do the mouse mean calculations and plot the mean points and lines
        
        %licks
        these_data_licks=lick_decision_times_pre(found_lick_decision_times_pre==1);
        these_mice_licks=these_mouse_nos_licks(found_lick_decision_times_pre==1);
         
        %peak
        these_data_peak=peak_decision_times_post_pref(found_peak_decision_times_post_pref==1);
        these_mice_peak=these_mouse_nos_peak(found_peak_decision_times_post_pref==1);
        
        %trough
        these_data_trough=trough_decision_times_post_pref(found_trough_decision_times_post_pref==1);
        these_mice_trough=these_mouse_nos_trough(found_trough_decision_times_post_pref==1);
        
        unique_mouse_nos=unique([these_mice_peak these_mice_trough these_mice_licks]);
        
        %Save data for correlation plots
        all_peak_dect=[all_peak_dect peak_decision_times_post_pref];
        all_peak_found_dect=[all_peak_found_dect found_peak_decision_times_post_pref];
        
        all_trough_dect=[all_trough_dect trough_decision_times_post_pref];
        all_trough_found_dect=[all_trough_found_dect found_trough_decision_times_post_pref];
        
        all_dect_groups=[all_dect_groups grNo*ones(1,length(found_peak_decision_times_post_pref))];
        
        if pacii==1
            all_lick_dect=[all_lick_dect lick_decision_times_pre];
            all_lick_found_dect=[all_lick_found_dect found_lick_decision_times_pre];
        end
        
         
        fprintf(1, ['Number of mice for ' group_legend{grNo} ' = ' num2str(length(unique_mouse_nos)) '\n\n']);
        fprintf(fileID, ['Number of mice for ' group_legend{grNo} ' = ' num2str(length(unique_mouse_nos)) '\n\n']);
        
        for msNo=unique_mouse_nos
            
            %peak
            this_mouse_mean_peak=mean(these_data_peak(these_mice_peak==msNo));
            
            glm_dect_mm.data(glm_ii_mm+1)=this_mouse_mean_peak;
            glm_dect_mm.lick_peak_trough(glm_ii_mm+1)=0;
            glm_dect_mm.pre_vs_post(glm_ii_mm+1)=1;
            glm_ii_mm=glm_ii_mm+1;
            
            %trough
            this_mouse_mean_trough=mean(these_data_trough(these_mice_trough==msNo));
            
            glm_dect_mm.data(glm_ii_mm+1)=this_mouse_mean_trough;
            glm_dect_mm.lick_peak_trough(glm_ii_mm+1)=1;
            glm_dect_mm.pre_vs_post(glm_ii_mm+1)=1;
            glm_ii_mm=glm_ii_mm+1;
            
            %licks
            this_mouse_mean_licks=mean(these_data_licks(these_mice_licks==msNo));
            
            glm_dect_mm.data(glm_ii_mm+1)=this_mouse_mean_licks;
            glm_dect_mm.lick_peak_trough(glm_ii_mm+1)=2;
            glm_dect_mm.pre_vs_post(glm_ii_mm+1)=1;
            glm_ii_mm=glm_ii_mm+1;
            
            
            plot([grNo grNo+4],[this_mouse_mean_peak this_mouse_mean_trough],...
                'o','MarkerSize',6,'MarkerFaceColor',[0.6 0.6 0.6],'MarkerEdgeColor',[0.6 0.6 0.6])
            
        end
        
        
        %Violin plot
        [mean_out, CIout]=drgViolinPoint(all_decision_dt_peak...
            ,edges,grNo,rand_offset,'k','k',3);
        
        %Violin plot
        [mean_out, CIout]=drgViolinPoint(all_decision_dt_trough...
            ,edges,grNo+4,rand_offset,'k','k',3);
        
%         %Violin plot
%         [mean_out, CIout]=drgViolinPoint(all_decision_dt_licks...
%             ,edges,grNo+8,rand_offset,'k','k',3);
        
        
        title(['Decision time post-lick ' bandwidth_names{pacii} ' prefrontal'])
        ylim([0 0.8])
        ylabel('dt')
        
        xticks([2 6])
        xticklabels({'Peak', 'Trough'})
        
    end
    
    %Perform the glm for decision time
    fprintf(1, ['glm for post-lick decision time per mouse per odor pair for prefrontal '  bandwidth_names{pacii} '\n'])
    fprintf(fileID, ['glm post-lick for decision time per mouse per odor pair for prefrontal '  bandwidth_names{pacii} '\n']);
    
    
    tbl = table(glm_dect.data',glm_dect.lick_peak_trough',glm_dect.group',...
        'VariableNames',{'detection_time','lick_peak_trough','genotype'});
    mdl = fitglm(tbl,'detection_time~genotype+lick_peak_trough+lick_peak_trough*genotype'...
        ,'CategoricalVars',[2,3])
    
    txt = evalc('mdl');
    txt=regexp(txt,'<strong>','split');
    txt=cell2mat(txt);
    txt=regexp(txt,'</strong>','split');
    txt=cell2mat(txt);
    
    fprintf(fileID,'%s\n', txt);
    
    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for post-lick detection time for prefrontal '  bandwidth_names{pacii} '\n'])
    fprintf(fileID, ['\n\nRanksum or t-test p values for post-lick detection time for prefrontal '  bandwidth_names{pacii} '\n']);
    
    
    [output_data] = drgMutiRanksumorTtest(input_dect_data, fileID);
    
    
    %Nested ANOVAN
    %https://www.mathworks.com/matlabcentral/answers/491365-within-between-subjects-in-anovan
    nesting=[0 0 0; ... % This line indicates that group factor is not nested in any other factor.
        0 0 0; ... % This line indicates that perCorr is not nested in any other factor.
        1 1 0];    % This line indicates that mouse_no (the last factor) is nested under group, perCorr and event
    % (the 1 in position 1 on the line indicates nesting under the first factor).
    figureNo = figureNo + 1;
    
    
    [p anovanTbl stats]=anovan(glm_dect.data,{glm_dect.lick_peak_trough glm_dect.group glm_dect.mouse_no},...
        'model','interaction',...
        'nested',nesting,...
        'varnames',{'lick_peak_trough','genotype','mouse_no'});
    
    fprintf(fileID, ['\n\nNested ANOVAN for post-lick decision time per mouse per odor pair for '  bandwidth_names{pacii} ' for prefrontal\n']);
    drgWriteANOVANtbl(anovanTbl,fileID);
    
    
    %Perform the glm for decision time to compare prefrontal vs. hippocampus
    fprintf(1, ['glm for post-lick decision time per mouse per odor pair for prefrontal and hippocampus '  bandwidth_names{pacii} '\n'])
    fprintf(fileID, ['glm post-lick for decision time per mouse per odor pair for prefrontal and hippocampus '  bandwidth_names{pacii} '\n']);
    
    
    tbl = table(glm_dect_all.data',glm_dect_all.lick_peak_trough',glm_dect_all.hipp_vs_pre',glm_dect_all.group',...
        'VariableNames',{'detection_time','lick_peak_trough','hippocampus_vs_prefrontal','genotype'});
    mdl = fitglm(tbl,'detection_time~genotype+lick_peak_trough+hippocampus_vs_prefrontal+lick_peak_trough*genotype*hippocampus_vs_prefrontal'...
        ,'CategoricalVars',[2,3,4])
    
    txt = evalc('mdl');
    txt=regexp(txt,'<strong>','split');
    txt=cell2mat(txt);
    txt=regexp(txt,'</strong>','split');
    txt=cell2mat(txt);
    
    fprintf(fileID,'%s\n', txt);
    
    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for post-lick detection time for hippocampus vs. prefrontal '  bandwidth_names{pacii} '\n'])
    fprintf(fileID, ['\n\nRanksum or t-test p values for post-lick detection time for hippocampus vs. prefrontal '  bandwidth_names{pacii} '\n']);
    
    
    [output_data] = drgMutiRanksumorTtest(input_dect_all_data, fileID);
    
    
    %Nested ANOVAN
    %https://www.mathworks.com/matlabcentral/answers/491365-within-between-subjects-in-anovan
    nesting=[0 0 0 0; ... % This line indicates that group factor is not nested in any other factor.
        0 0 0 0; ... % This line indicates that perCorr is not nested in any other factor.
        0 0 0 0; ... % This line indicates that perCorr is not nested in any other factor.
        1 1  1 0];    % This line indicates that mouse_no (the last factor) is nested under group, perCorr and event
    % (the 1 in position 1 on the line indicates nesting under the first factor).
    figureNo = figureNo + 1;
    
    
    [p anovanTbl stats]=anovan(glm_dect_all.data,{glm_dect_all.lick_peak_trough glm_dect_all.group glm_dect_all.hipp_vs_pre glm_dect_all.mouse_no},...
        'model','interaction',...
        'nested',nesting,...
        'varnames',{'lick_peak_trough','genotype','hippocampus_vs_prefrontal','mouse_no'});
    
    fprintf(fileID, ['\n\nNested ANOVAN for post-lick decision time per mouse per odor pair for '  bandwidth_names{pacii} ' for hippocampus vs. prefrontal\n']);
    drgWriteANOVANtbl(anovanTbl,fileID);
    
    %Perform the glm for decision time per mouse
    fprintf(1, ['glm for decision time per mouse for prefrontal\n'])
    %     fprintf(fileID, ['glm for decision time per mouse for prefrontal\n']);
    
    tbl = table(glm_dect_mm.data',glm_dect_mm.lick_peak_trough',glm_dect_mm.pre_vs_post',...
        'VariableNames',{'detection_time','lick_peak_trough','pre_vs_post'});
    mdl = fitglm(tbl,'detection_time~lick_peak_trough+pre_vs_post+lick_peak_trough*pre_vs_post'...
        ,'CategoricalVars',[2,3])
    
  %Plot the correlations for the peak
    figureNo = figureNo + 1;
    try
        close(figureNo)
    catch
    end
    hFig=figure(figureNo);
    
    set(hFig, 'units','normalized','position',[.1 .5 .25 .4])
    
    hold on
    
    ax=gca;ax.LineWidth=3;
    
    for grNo=1:2
        these_included=all_peak_found_dect&all_lick_found_dect&(all_dect_groups==grNo);
        switch grNo
            case 1
                plot(all_peak_dect(these_included),all_lick_dect(these_included),'o','MarkerEdgeColor','g','MarkerFaceColor','g')
            case 2
                plot(all_peak_dect(these_included),all_lick_dect(these_included),'o','MarkerEdgeColor','b','MarkerFaceColor','b')
            case 3
                plot(all_peak_dect(these_included),all_lick_dect(these_included),'o','MarkerEdgeColor','y','MarkerFaceColor','y')
        end
        
    end
    
    these_included=all_peak_found_dect&all_lick_found_dect;
    y=all_lick_dect(these_included)';
    x=all_peak_dect(these_included)';
    
    
    [rho,pval]=corr(x,y);
    
    %Fit a line
    c = polyfit(x,y,1);
    % Display evaluated equation y = m*x + b
    disp(['Equation is y = ' num2str(c(1)) '*x + ' num2str(c(2))])
    % Evaluate fit equation using polyval
    y_est = polyval(c,x);
    % Add trend line to plot
    plot(x,y_est,'k-','LineWidth',2)
    
    xlim([0 0.8])
    ylim([0 0.8])
    
    fprintf(1, ['\nrho = %d, p value = %d for lick discrimination time vs prefrontal peak discrimination time for '  bandwidth_names{pacii} '\n\n'],rho,pval)
    
    ylabel('Discrimintion time for licks')
    xlabel('Discrimination time for peak pre-lick ztPRP')
    title(['Discrimination time for licks vs peak post-lick ztPRP '  bandwidth_names{pacii} ' prefrontal' ])
    
    %Plot the correlations for the trough
    figureNo = figureNo + 1;
    try
        close(figureNo)
    catch
    end
    hFig=figure(figureNo);
    
    set(hFig, 'units','normalized','position',[.1 .5 .25 .4])
    
    hold on
    
    ax=gca;ax.LineWidth=3;
    
    for grNo=1:2
        these_included=all_trough_found_dect&all_lick_found_dect&(all_dect_groups==grNo);
        switch grNo
            case 1
                plot(all_trough_dect(these_included),all_lick_dect(these_included),'o','MarkerEdgeColor','g','MarkerFaceColor','g')
            case 2
                plot(all_trough_dect(these_included),all_lick_dect(these_included),'o','MarkerEdgeColor','b','MarkerFaceColor','b')
            case 3
                plot(all_trough_dect(these_included),all_lick_dect(these_included),'o','MarkerEdgeColor','y','MarkerFaceColor','y')
        end
        
    end
    
    these_included=all_trough_found_dect&all_lick_found_dect;
    y=all_lick_dect(these_included)';
    x=all_trough_dect(these_included)';
    
    
    [rho,pval]=corr(x,y);
    
    %Fit a line
    c = polyfit(x,y,1);
    % Display evaluated equation y = m*x + b
    disp(['Equation is y = ' num2str(c(1)) '*x + ' num2str(c(2))])
    % Evaluate fit equation using polyval
    y_est = polyval(c,x);
    % Add trend line to plot
    plot(x,y_est,'k-','LineWidth',2)
    
    xlim([0 0.8])
    ylim([0 0.8])
    
    fprintf(1, ['\nrho = %d, p value = %d for lick discrimination time vs prefrontal trough ztPRP discrimination time for '  bandwidth_names{pacii} '\n\n'],rho,pval)
    
    ylabel('Discrimintion time for licks')
    xlabel('Discrimination time for trough pre-lick ztPRP')
    title(['Discrimination time for licks vs trough post-lick ztPRP '  bandwidth_names{pacii} ' prefrontal' ])
end



fclose(fileID);


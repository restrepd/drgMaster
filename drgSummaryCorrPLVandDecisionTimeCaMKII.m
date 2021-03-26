function drgSummaryCorrPLVandDecisionTimeCaMKII
%Analyzes the linear discriminant analysis performed by drgLFPDiscriminantBatch
%Takes as a in input the 'drgbDiscPar' file listing 'Discriminant_*.mat' output files
%from drgLFPDiscriminantBatch
%
%Performs summary analises for LDA and  PAC analysis performed with case 19
%of drgAnalysisBatchLFP

warning('off')

close all
clear all


bandwidth_names{1}='Theta';
bandwidth_names{2}='Beta';
bandwidth_names{3}='Low gamma';
bandwidth_names{4}='High gamma';

prof_naive_leg{1}='Proficient';
prof_naive_leg{2}='Naive';

group_legend{1}='WT';
group_legend{2}='Het';
group_legend{3}='KO';

evTypeLabels{1}='S+';
evTypeLabels{2}='S-';

peak_label{1}='Trough';
peak_label{2}='Peak';

%Location of files
% hippPathName='E:\CaMKIIpaper\datos sumarry\coherence\';
PLV_PathName='/Users/restrepd/Documents/Projects/CaMKII_analysis/PLV/';

%Files
PLV_FileName{1}='CaMKIIPLVCamKIAPEB0301202_out.mat';
PLV_FileName{2}='CaMKIIEBAPPLV02282021_out.mat';
PLV_FileName{3}='CaMKIIPAEAPLV03062021_out.mat';
PLV_FileName{4}='CaMKIIEAPAPLV03112021_out.mat';
PLV_FileName{5}='CaMKIIpz1EAPAPLV03052021_out.mat';
PLV_FileName{6}='CaMKIIpz1PAEAPLV03042021_out.mat';
PLV_FileName{7}='CaMKIIpzz1EAPAPLV02262021_out.mat';
PLV_FileName{8}='CaMKIIpzz1propylacecPLV02222021_out.mat';

hippPathName='/Users/restrepd/Documents/Projects/CaMKII_analysis/Discriminant/';

%Files for decision times

%Hippocampus

%pzz1PAEA
hippFileName{8}='pcorr_Discriminant_CaMKIIpzz1paea_disc_01302021_hippo2.mat';
file_legend{8}='pzz1PAEA';

%pzz1EAPA
hippFileName{7}='pcorr_Discriminant_CaMKIIpzz1ethylace_disc_PRPhipp2_01222020.mat';
file_legend{7}='pzz1EAPA';

%pz1PAEA
hippFileName{6}='pcorr_Discriminant_CaMKIIpz1paea_disc_02042021_hipp2.mat';
file_legend{6}='pz1PAEA';

%pz1EAPA
hippFileName{5}='pcorr_Discriminant_CaMKIIPZ1EAPA_disc_PRPhippo2_0212021.mat';
file_legend{5}='pz1EAPA';

%PAEA
hippFileName{4}='pcorr_Discriminant_CaMKIIPAEA_disc_PRPhipp2_01222020.mat';
file_legend{4}='PAEA';

%EBAP
hippFileName{2}='pcorr_Discriminant_CaMKIIEBAP_disc_PRPhipp2_01222020.mat';
file_legend{2}='EBAP';

%EAPA
hippFileName{3}='pcorr_Discriminant_CaMKIIEAPA_disc_PRPhippo2_02082021.mat';
file_legend{3}='EAPA';

%APEB aceto
hippFileName{1}='pcorr_Discriminant_CaMKIIaceto_disc_PRPhippo2_0272021.mat';
file_legend{1}='APEB';


%Load data
all_PLV_files=[];
all_beh_files=[];
figNo=0;

%Now
for bwii=1:4    %for amplitude bandwidths (beta, low gamma, high gamma)
    
    all_ddPLV=[];
    all_disc_times=[];
    
    id_ii=0;
    input_data=[];
    
    %Plot the average
    figNo = figNo +1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    hold on
    
    for ii=1:length(PLV_FileName)
        load([PLV_PathName PLV_FileName{ii}])
%         all_PLV_files(ii).delta_PLV_odor_minus_ref_per_mouse=delta_PLV_odor_minus_ref_per_mouse;
%         all_PLV_files(ii).delta_phase_odor_per_mouse=delta_phase_odor_per_mouse;
%         all_PLV_files(ii).group_no_per_mouse=group_no_per_mouse;
        PLV_group_no_per_mouse=group_no_per_mouse;
        load([hippPathName hippFileName{ii}])
%         all_beh_files(ii).mean_per_corr_per_mouse_prof=mean_per_corr_per_mouse_prof;
%         all_beh_files(ii).group_no_per_mouse=group_no_per_mouse;
        
        per_ii=1;
        
        
        
        for grNo=1:3
            pfft=1;
            these_mice=disc_time.PACii(1).pcorr(per_ii).group(grNo).mouseNos;
            these_ddPLV=[];
            these_disc_times=[];
            for mouseNo=1:length(PLV_group_no_per_mouse)
                %PLV inlcudes all mice, but decision time may not
                if PLV_group_no_per_mouse(mouseNo)==grNo
                    if sum(these_mice==mouseNo)>0
                        this_ddPLV=delta_PLV_odor_minus_ref_per_mouse(mouseNo,bwii,per_ii,1)-delta_PLV_odor_minus_ref_per_mouse(mouseNo,bwii,per_ii,2);
                        all_ddPLV=[all_ddPLV this_ddPLV'];
                        these_ddPLV=[these_ddPLV this_ddPLV];
                        ii_mouse=find(these_mice==mouseNo);
                        this_disc_time=disc_time.PACii(1).licks.pcorr(per_ii).group(grNo).data(ii_mouse);
                        all_disc_times=[all_disc_times this_disc_time];
                        these_disc_times=[these_disc_times this_disc_time];
                    end
                end
            end
            
            switch grNo
                case 1
                    plot(these_ddPLV,these_disc_times,'o','MarkerEdgeColor','k','MarkerFaceColor','g')
                case 2
                    plot(these_ddPLV,these_disc_times,'o','MarkerEdgeColor','k','MarkerFaceColor','b')
                case 3
                    plot(these_ddPLV,these_disc_times,'o','MarkerEdgeColor','k','MarkerFaceColor','y')
            end
            
            
        end
        
    end
    [rho,pval(bwii)]=corr(all_ddPLV(~isnan(all_disc_times))',all_disc_times(~isnan(all_disc_times))');
    
    fprintf(1, ['rho = %d, p value = %d for '  bandwidth_names{bwii} '\n'],rho,pval(bwii))
    
    xlabel('ddPLV')
    ylabel('Percent correct')
    title(bandwidth_names{bwii} )
end
 
 fprintf(1,'\n')
 
for bwii=1:4    %for amplitude bandwidths (beta, low gamma, high gamma)
    
    all_ddphase=[];
    all_disc_times=[];
    
    id_ii=0;
    input_data=[];
    
    %Plot the average
    figNo = figNo +1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    hold on
    
    for ii=1:length(PLV_FileName)
        load([PLV_PathName PLV_FileName{ii}])
%         all_PLV_files(ii).delta_PLV_odor_minus_ref_per_mouse=delta_PLV_odor_minus_ref_per_mouse;
%         all_PLV_files(ii).delta_phase_odor_per_mouse=delta_phase_odor_per_mouse;
%         all_PLV_files(ii).group_no_per_mouse=group_no_per_mouse;
        PLV_group_no_per_mouse=group_no_per_mouse;
        load([hippPathName hippFileName{ii}])
%         all_beh_files(ii).mean_per_corr_per_mouse_prof=mean_per_corr_per_mouse_prof;
%         all_beh_files(ii).group_no_per_mouse=group_no_per_mouse;
        
        per_ii=1;
        
        
        for grNo=1:3
            pfft=1;
            
               these_mice=disc_time.PACii(1).pcorr(per_ii).group(grNo).mouseNos;
            these_ddphase=[];
            these_disc_times=[];
            for mouseNo=1:length(PLV_group_no_per_mouse)
                %PLV inlcudes all mice, but decision time may not
                if PLV_group_no_per_mouse(mouseNo)==grNo
                    if sum(these_mice==mouseNo)>0
                        this_ddphase=delta_phase_odor_per_mouse(mouseNo,bwii,per_ii,1)-delta_PLV_odor_minus_ref_per_mouse(mouseNo,bwii,per_ii,2);
                        all_ddphase=[all_ddphase this_ddphase'];
                        these_ddphase=[these_ddphase this_ddphase];
                        ii_mouse=find(these_mice==mouseNo);
                        this_disc_time=disc_time.PACii(1).licks.pcorr(per_ii).group(grNo).data(ii_mouse);
                        all_disc_times=[all_disc_times this_disc_time];
                        these_disc_times=[these_disc_times this_disc_time];
                    end
                end
            end
            
            switch grNo
                case 1
                    plot(these_ddphase,these_disc_times,'o','MarkerEdgeColor','k','MarkerFaceColor','g')
                case 2
                    plot(these_ddphase,these_disc_times,'o','MarkerEdgeColor','k','MarkerFaceColor','b')
                case 3
                    plot(these_ddphase,these_disc_times,'o','MarkerEdgeColor','k','MarkerFaceColor','y')
            end
            
            
        end
        
    end
    [rho,pval]=corr(all_ddphase(~isnan(all_disc_times))',all_disc_times(~isnan(all_disc_times))');
    
    fprintf(1, ['rho = %d, p value = %d for '  bandwidth_names{bwii} '\n'],rho,pval)
    
    xlabel('ddphase')
    ylabel('Percent correct')
    title(bandwidth_names{bwii} )
end


%Load data



%Now do p value
for bwii=1:4    %for amplitude bandwidths (beta, low gamma, high gamma)
    
    all_ddPLV=[];
    all_pval_outs=[];
    
    id_ii=0;
    input_data=[];
    
    %Plot the average
    figNo = figNo +1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    hold on
    
    for ii=1:length(PLV_FileName)
        load([PLV_PathName PLV_FileName{ii}])
%         all_PLV_files(ii).delta_PLV_odor_minus_ref_per_mouse=delta_PLV_odor_minus_ref_per_mouse;
%         all_PLV_files(ii).delta_phase_odor_per_mouse=delta_phase_odor_per_mouse;
%         all_PLV_files(ii).group_no_per_mouse=group_no_per_mouse;
        PLV_group_no_per_mouse=group_no_per_mouse;
        load([hippPathName hippFileName{ii}])
%         all_beh_files(ii).mean_per_corr_per_mouse_prof=mean_per_corr_per_mouse_prof;
%         all_beh_files(ii).group_no_per_mouse=group_no_per_mouse;
        
        per_ii=1;
        
        
        
        for grNo=1:3
            pfft=1;
            these_mice=pval_out.PACii(1).pcorr(per_ii).group(grNo).mouseNos;
            these_ddPLV=[];
            these_pval_outs=[];
            for mouseNo=1:length(PLV_group_no_per_mouse)
                %PLV inlcudes all mice, but decision time may not
                if PLV_group_no_per_mouse(mouseNo)==grNo
                    if sum(these_mice==mouseNo)>0
                        this_ddPLV=delta_PLV_odor_minus_ref_per_mouse(mouseNo,bwii,per_ii,1)-delta_PLV_odor_minus_ref_per_mouse(mouseNo,bwii,per_ii,2);
                        all_ddPLV=[all_ddPLV this_ddPLV'];
                        these_ddPLV=[these_ddPLV this_ddPLV];
                        ii_mouse=find(these_mice==mouseNo);
                        this_pval_out=pval_out.PACii(1).lick.pcorr(per_ii).group(grNo).odor_data(ii_mouse);
                        all_pval_outs=[all_pval_outs this_pval_out];
                        these_pval_outs=[these_pval_outs this_pval_out];
                    end
                end
            end
            
            switch grNo
                case 1
                    plot(these_ddPLV,these_pval_outs,'o','MarkerEdgeColor','k','MarkerFaceColor','g')
                case 2
                    plot(these_ddPLV,these_pval_outs,'o','MarkerEdgeColor','k','MarkerFaceColor','b')
                case 3
                    plot(these_ddPLV,these_pval_outs,'o','MarkerEdgeColor','k','MarkerFaceColor','y')
            end
            
            
        end
        
    end
    [rho,pval(bwii)]=corr(all_ddPLV(~isnan(all_pval_outs))',all_pval_outs(~isnan(all_pval_outs))');
    
    fprintf(1, ['rho = %d, p value = %d for '  bandwidth_names{bwii} '\n'],rho,pval(bwii))
    
    xlabel('ddPLV')
    ylabel('log(pval)')
    title(bandwidth_names{bwii} )
end
 
 fprintf(1,'\n')
 
for bwii=1:4    %for amplitude bandwidths (beta, low gamma, high gamma)
    
    all_ddphase=[];
    all_pval_outs=[];
    
    id_ii=0;
    input_data=[];
    
    %Plot the average
    figNo = figNo +1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    hold on
    
    for ii=1:length(PLV_FileName)
        load([PLV_PathName PLV_FileName{ii}])
%         all_PLV_files(ii).delta_PLV_odor_minus_ref_per_mouse=delta_PLV_odor_minus_ref_per_mouse;
%         all_PLV_files(ii).delta_phase_odor_per_mouse=delta_phase_odor_per_mouse;
%         all_PLV_files(ii).group_no_per_mouse=group_no_per_mouse;
        PLV_group_no_per_mouse=group_no_per_mouse;
        load([hippPathName hippFileName{ii}])
%         all_beh_files(ii).mean_per_corr_per_mouse_prof=mean_per_corr_per_mouse_prof;
%         all_beh_files(ii).group_no_per_mouse=group_no_per_mouse;
        
        per_ii=1;
        
        
        for grNo=1:3
            pfft=1;
            
               these_mice=pval_out.PACii(1).pcorr(per_ii).group(grNo).mouseNos;
            these_ddphase=[];
            these_pval_outs=[];
            for mouseNo=1:length(PLV_group_no_per_mouse)
                %PLV inlcudes all mice, but decision time may not
                if PLV_group_no_per_mouse(mouseNo)==grNo
                    if sum(these_mice==mouseNo)>0
                        this_ddphase=delta_phase_odor_per_mouse(mouseNo,bwii,per_ii,1)-delta_PLV_odor_minus_ref_per_mouse(mouseNo,bwii,per_ii,2);
                        all_ddphase=[all_ddphase this_ddphase'];
                        these_ddphase=[these_ddphase this_ddphase];
                        ii_mouse=find(these_mice==mouseNo);
                        this_pval_out=pval_out.PACii(1).lick.pcorr(per_ii).group(grNo).odor_data(ii_mouse);
                        all_pval_outs=[all_pval_outs this_pval_out];
                        these_pval_outs=[these_pval_outs this_pval_out];
                    end
                end
            end
            
            switch grNo
                case 1
                    plot(these_ddphase,these_pval_outs,'o','MarkerEdgeColor','k','MarkerFaceColor','g')
                case 2
                    plot(these_ddphase,these_pval_outs,'o','MarkerEdgeColor','k','MarkerFaceColor','b')
                case 3
                    plot(these_ddphase,these_pval_outs,'o','MarkerEdgeColor','k','MarkerFaceColor','y')
            end
            
            
        end
        
    end
    [rho,pval]=corr(all_ddphase(~isnan(all_pval_outs))',all_pval_outs(~isnan(all_pval_outs))');
    
    fprintf(1, ['rho = %d, p value = %d for '  bandwidth_names{bwii} '\n'],rho,pval)
    
    xlabel('ddphase')
    ylabel('log(pval)')
    title(bandwidth_names{bwii} )
end
pffft=1;
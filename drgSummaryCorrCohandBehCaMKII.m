function drgSummaryCorrCohandBehCaMKII
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

behPathName='/Users/restrepd/Documents/Projects/CaMKII_analysis/Behavior/';

%Behavioral files
behFileName{1}='CaMKIICamKIAPEB_beh_03142021beh_80_out.mat';
behFileName{2}='CaMKIIEAPAbeh03142021beh_80_out.mat';
behFileName{3}='CaMKIIEBAP_beh_03132021beh_80_out.mat';
behFileName{4}='CaMKIIPAEA_beh_03132021beh_80_out.mat';
behFileName{5}='CaMKIIpz1eapa_beh_03132021beh_80_out.mat';
behFileName{6}='CaMKIIpz1paea_beh_03132021beh_80_out.mat';
behFileName{7}='CaMKIIpzz1ethylace_beh_03132021beh_80_out.mat';
behFileName{8}='CaMKIIpzz1paea_beh_04022021beh_80_out.mat';

%Coherence files
coh_PathName='/Users/restrepd/Documents/Projects/CaMKII_analysis/Coherence new/';

%Files
coh_FileName{1}='CaMKIIacetocohe02012021_out80.mat';
coh_FileName{2}='CaMKIIEAPAcohe2262021_out80.mat';
coh_FileName{3}='CaMKIIethylbenacetocohe2262021_out80.mat';
coh_FileName{4}='CaMKIIPAEAcohe02082021_out80.mat';
coh_FileName{5}='CaMKIIpz1EAcohe02142021_out80.mat';
coh_FileName{6}='CaMKIIPZ1PAEAcohe202102021_out80.mat';
coh_FileName{7}='CaMKIIpzz1EAPAcohe02112021_out80.mat';
coh_FileName{8}='CaMKIIpzz1propylacecohe02092021_out80.mat';

%Correlation for coherence

%Load data
all_coh_files=[];
all_beh_files=[];
figNo=0;

%Now
for bwii=1:4    %for amplitude bandwidths (beta, low gamma, high gamma)
    
    all_coh_auROC=[];
    all_pCorr=[];
    
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
     
    for ii=1:length(coh_FileName)
        load([coh_PathName coh_FileName{ii}])
%         all_PLV_files(ii).delta_PLV_odor_minus_ref_per_mouse=delta_PLV_odor_minus_ref_per_mouse;
%         all_PLV_files(ii).delta_phase_odor_per_mouse=delta_phase_odor_per_mouse;
%         all_PLV_files(ii).group_no_per_mouse=group_no_per_mouse;
        
        load([behPathName behFileName{ii}])
%         all_beh_files(ii).mean_per_corr_per_mouse_prof=mean_per_corr_per_mouse_prof;
%         all_beh_files(ii).group_no_per_mouse=group_no_per_mouse;
        
        per_ii=1;
        
        
        for grNo=1:3
            pfft=1;
            
            these_coh_auROC=[];
            these_coh_auROC=coh_auROC_per_mouse(coh_group_no_per_mouse==grNo,bwii,per_ii);
            all_coh_auROC=[all_coh_auROC these_coh_auROC'];
            these_pCorr=mean_per_corr_per_mouse_prof(group_no_per_mouse==grNo);
            all_pCorr=[all_pCorr these_pCorr];
            switch grNo
                case 1
                    plot(these_coh_auROC,these_pCorr,'o','MarkerEdgeColor','k','MarkerFaceColor','g')
                case 2
                    plot(these_coh_auROC,these_pCorr,'o','MarkerEdgeColor','k','MarkerFaceColor','b')
                case 3
                    plot(these_coh_auROC,these_pCorr,'o','MarkerEdgeColor','k','MarkerFaceColor','y')
            end
            
            
        end
        
    end
    [rho,pval(bwii)]=corr(all_coh_auROC(~isnan(all_pCorr))',all_pCorr(~isnan(all_pCorr))');
    
    fprintf(1, ['rho = %d, p value = %d for '  bandwidth_names{bwii} '\n'],rho,pval(bwii))
    
    xlabel('auROC')
    ylabel('Percent correct')
    title(['Coherence auROC vs percent correct' bandwidth_names{bwii}])
end
 
 fprintf(1,'\n')

%Correlation for PLV
%Load data
all_PLV_files=[];
all_beh_files=[];

%Now
for bwii=1:4    %for amplitude bandwidths (beta, low gamma, high gamma)
    
    all_ddPLV=[];
    all_pCorr=[];
    
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
        load([behPathName behFileName{ii}])
%         all_beh_files(ii).mean_per_corr_per_mouse_prof=mean_per_corr_per_mouse_prof;
%         all_beh_files(ii).group_no_per_mouse=group_no_per_mouse;
        
        per_ii=1;
        
        
        for grNo=1:3
            pfft=1;
            
            these_ddPLV=delta_PLV_odor_minus_ref_per_mouse(PLV_group_no_per_mouse==grNo,bwii,per_ii,1)-delta_PLV_odor_minus_ref_per_mouse(group_no_per_mouse==grNo,bwii,per_ii,2);
            all_ddPLV=[all_ddPLV these_ddPLV'];
            these_pCorr=mean_per_corr_per_mouse_prof(group_no_per_mouse==grNo);
            all_pCorr=[all_pCorr these_pCorr];
            switch grNo
                case 1
                    plot(these_ddPLV,these_pCorr,'o','MarkerEdgeColor','k','MarkerFaceColor','g')
                case 2
                    plot(these_ddPLV,these_pCorr,'o','MarkerEdgeColor','k','MarkerFaceColor','b')
                case 3
                    plot(these_ddPLV,these_pCorr,'o','MarkerEdgeColor','k','MarkerFaceColor','y')
            end
            
            
        end
        
    end
    [rho,pval(bwii)]=corr(all_ddPLV(~isnan(all_pCorr))',all_pCorr(~isnan(all_pCorr))');
    
    fprintf(1, ['rho = %d, p value = %d for '  bandwidth_names{bwii} '\n'],rho,pval(bwii))
    
    xlabel('ddPLV')
    ylabel('Percent correct')
    title(bandwidth_names{bwii} )
end
 
 fprintf(1,'\n')
 
for bwii=1:4    %for amplitude bandwidths (beta, low gamma, high gamma)
    
    all_ddphase=[];
    all_pCorr=[];
    
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
        load([behPathName behFileName{ii}])
%         all_beh_files(ii).mean_per_corr_per_mouse_prof=mean_per_corr_per_mouse_prof;
%         all_beh_files(ii).group_no_per_mouse=group_no_per_mouse;
        
        per_ii=1;
        
        
        for grNo=1:3
            pfft=1;
            
            these_ddphase=delta_phase_odor_per_mouse(PLV_group_no_per_mouse==grNo,bwii,per_ii,1);
            all_ddphase=[all_ddphase these_ddphase'];
            these_pCorr=mean_per_corr_per_mouse_prof(group_no_per_mouse==grNo);
            all_pCorr=[all_pCorr these_pCorr];
            switch grNo
                case 1
                    plot(these_ddphase,these_pCorr,'o','MarkerEdgeColor','k','MarkerFaceColor','g')
                case 2
                    plot(these_ddphase,these_pCorr,'o','MarkerEdgeColor','k','MarkerFaceColor','b')
                case 3
                    plot(these_ddphase,these_pCorr,'o','MarkerEdgeColor','k','MarkerFaceColor','y')
            end
            
            
        end
        
    end
    [rho,pval]=corr(all_ddphase(~isnan(all_pCorr))',all_pCorr(~isnan(all_pCorr))');
    
    fprintf(1, ['rho = %d, p value = %d for '  bandwidth_names{bwii} '\n'],rho,pval)
    
    xlabel('ddphase')
    ylabel('Percent correct')
    title(bandwidth_names{bwii} )
end


for bwii=1:4    %for amplitude bandwidths (beta, low gamma, high gamma)
    
    all_ddphase=[];
    all_ddPLV=[];
    
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
        load([behPathName behFileName{ii}])
%         all_beh_files(ii).mean_per_corr_per_mouse_prof=mean_per_corr_per_mouse_prof;
%         all_beh_files(ii).group_no_per_mouse=group_no_per_mouse;
        
        per_ii=1;
        
        
        for grNo=1:3
            pfft=1;
            
             these_ddPLV=delta_PLV_odor_minus_ref_per_mouse(PLV_group_no_per_mouse==grNo,bwii,per_ii,1)-delta_PLV_odor_minus_ref_per_mouse(group_no_per_mouse==grNo,bwii,per_ii,2);
            all_ddPLV=[all_ddPLV these_ddPLV'];
            these_ddphase=delta_phase_odor_per_mouse(PLV_group_no_per_mouse==grNo,bwii,per_ii,1);
            all_ddphase=[all_ddphase these_ddphase'];
           
            switch grNo
                case 1
                    plot(these_ddphase,these_ddPLV,'o','MarkerEdgeColor','k','MarkerFaceColor','g')
                case 2
                    plot(these_ddphase,these_ddPLV,'o','MarkerEdgeColor','k','MarkerFaceColor','b')
                case 3
                    plot(these_ddphase,these_ddPLV,'o','MarkerEdgeColor','k','MarkerFaceColor','y')
            end
            
            
        end
        
    end
    [rho,pval]=corr(all_ddphase',all_ddPLV');
    
    fprintf(1, ['rho = %d, p value = %d for '  bandwidth_names{bwii} '\n'],rho,pval)
    
    xlabel('ddphase')
    ylabel('ddPLV')
    title(bandwidth_names{bwii} )
end


pffft=1;
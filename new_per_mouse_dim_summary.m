function new_per_mouse_dim_summary


warning('off')

close all
clear all
 
no_odor_pairs=6;

odor_pair_no=0;



no_odor_pairs=6;
 
odor_pair_no=0;

%Daniel's
outFileName='Discriminant_spmc_discriminantolfac_all_PCA_LFP_aceto42919.mat';
try
    outPathName='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 5 LDA/APEBDR/';
    load([outPathName outFileName])
catch
    outPathName='C:\Users\Justin Losacco\Dropbox\Restrepo_Lab\Papers\JL paper\Identity paper\SPM manuscript  Justin-Daniel\Fig 5 LDA\APEBDR\';
    load([outPathName outFileName])
end
odor_pair_no=odor_pair_no+1;
handles_dim_per_odor_pair(odor_pair_no).handles_out=handles_out;
odor_pair_label{1}='APEBloc1';
exp(1)=1;
 
outFileName='Discriminant_spmc_discriminantolfac_all_PCA_LFP_EAPA42919.mat';
try
    outPathName='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 5 LDA/EAPADR/';
    load([outPathName outFileName])
catch
    outPathName='C:\Users\Justin Losacco\Dropbox\Restrepo_Lab\Papers\JL paper\Identity paper\SPM manuscript  Justin-Daniel\Fig 5 LDA\EAPADR\';
    load([outPathName outFileName])
end
odor_pair_no=odor_pair_no+1;
handles_dim_per_odor_pair(odor_pair_no).handles_out=handles_out;
odor_pair_label{2}='EAPAloc1';
exp(2)=1;

outFileName='Discriminant_spmc_discriminantolfac_all_PCA_LFP_iso42919.mat';
try
    outPathName='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 5 LDA/IAMODR/';
    load([outPathName outFileName])
catch
    outPathName='C:\Users\Justin Losacco\Dropbox\Restrepo_Lab\Papers\JL paper\Identity paper\SPM manuscript  Justin-Daniel\Fig 5 LDA\IAMODR\';
    load([outPathName outFileName])
end
odor_pair_no=odor_pair_no+1;
handles_dim_per_odor_pair(odor_pair_no).handles_out=handles_out;
odor_pair_label{3}='IAMOloc1';
exp(3)=1;

%Justin's
outFileName='Discriminant_spm_discriminant_LFP_wavephase_05172019_IAAP.mat';
try
    outPathName='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 5 LDA/IAAPJL/';
    load([outPathName outFileName])
catch
    outPathName='C:\Users\Justin Losacco\Dropbox\Restrepo_Lab\Papers\JL paper\Identity paper\SPM manuscript  Justin-Daniel\Fig 5 LDA\IAAPJL\';
    load([outPathName outFileName])
end
odor_pair_no=odor_pair_no+1;
handles_dim_per_odor_pair(odor_pair_no).handles_out=handles_out;
odor_pair_label{4}='IAAPloc2';
exp(4)=2;

outFileName='Discriminant_spm_discriminant_LFP_wavephase_04252019_EAPA.mat';
try
    outPathName='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 5 LDA/EAPAJL/';
    load([outPathName outFileName])
catch
    outPathName='C:\Users\Justin Losacco\Dropbox\Restrepo_Lab\Papers\JL paper\Identity paper\SPM manuscript  Justin-Daniel\Fig 5 LDA\EAPAJL\';
    load([outPathName outFileName])
end
odor_pair_no=odor_pair_no+1;
handles_dim_per_odor_pair(odor_pair_no).handles_out=handles_out;
odor_pair_label{5}='EAPAloc2';
exp(5)=2;

outFileName='Discriminant_spm_discriminant_LFP_wavephase_04252019_IsoAA_mo.mat';
try
    outPathName='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 5 LDA/IAMOJL/';
    load([outPathName outFileName])
catch
    outPathName='C:\Users\Justin Losacco\Dropbox\Restrepo_Lab\Papers\JL paper\Identity paper\SPM manuscript  Justin-Daniel\Fig 5 LDA\IAMOJL\';
    load([outPathName outFileName])
end
odor_pair_no=odor_pair_no+1;
handles_dim_per_odor_pair(odor_pair_no).handles_out=handles_out;
odor_pair_label{6}='IAMOloc2';
exp(6)=2;

PACnames{1}='Beta';
PACnames{2}='Low gamma';
PACnames{3}='High gamma';

window_label{1}='Pre-odor';
window_label{2}='Odor';

% 
% 
% 
% %Daniel's
% dimFileName='Discriminant_2spm_discriminantolfac_per_mice_aceto060519.mat';
% try
%     dimPathName='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 5 Decision time-LDA all mouse/APEB DR AM/';
%     load([dimPathName dimFileName])
% catch
%     dimPathName='C:\Users\Justin Losacco\Dropbox\Restrepo_Lab\Papers\JL paper\Identity paper\SPM manuscript  Justin-Daniel\Fig 5 Decision time-LDA all mouse\APEB DR AM\';
%     load([dimPathName dimFileName])
% end
% odor_pair_no=odor_pair_no+1;
% handles_dim_per_odor_pair(odor_pair_no).handles_out=handles_out;
% odor_pair_label{1}='APEBloc1';
% exp(1)=1;
% 
% dimFileName='Discriminant_2spmc_discriminantolfac_per_mice_EAPA06062019.mat';
% try
%     dimPathName='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 5 Decision time-LDA all mouse/EAPA DR AM/';
%     load([dimPathName dimFileName])
% catch
%     dimPathName='C:\Users\Justin Losacco\Dropbox\Restrepo_Lab\Papers\JL paper\Identity paper\SPM manuscript  Justin-Daniel\Fig 5 Decision time-LDA all mouse\EAPA DR AM\';
%     load([dimPathName dimFileName])
% end
% odor_pair_no=odor_pair_no+1;
% handles_dim_per_odor_pair(odor_pair_no).handles_out=handles_out;
% odor_pair_label{2}='EAPAloc1';
% exp(2)=1;
% 
% dimFileName='Discriminant_2spm_discriminantolfac_per_mice_iso06082019.mat';
% try
%     dimPathName='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 5 Decision time-LDA all mouse/IAMO DR AM/';
%     load([dimPathName dimFileName])
% catch
%     dimPathName='C:\Users\Justin Losacco\Dropbox\Restrepo_Lab\Papers\JL paper\Identity paper\SPM manuscript  Justin-Daniel\Fig 5 Decision time-LDA all mouse\IAMO DR AM\';
%     load([dimPathName dimFileName])
% end
% odor_pair_no=odor_pair_no+1;
% handles_dim_per_odor_pair(odor_pair_no).handles_out=handles_out;
% odor_pair_label{3}='IAMOloc1';
% exp(3)=1;
% 
% %Justin's
% dimFileName='Discriminant_2spm_discriminant_LFP_wp_per_mouse_05312019_IAAP.mat';
% try
%     dimPathName='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 5 Decision time-LDA all mouse/IAAP JL AM/';
%     load([dimPathName dimFileName])
% catch
%     dimPathName='C:\Users\Justin Losacco\Dropbox\Restrepo_Lab\Papers\JL paper\Identity paper\SPM manuscript  Justin-Daniel\Fig 5 Decision time-LDA all mouse\IAAP JL AM\';
%     load([dimPathName dimFileName])
% end
% odor_pair_no=odor_pair_no+1;
% handles_dim_per_odor_pair(odor_pair_no).handles_out=handles_out;
% odor_pair_label{4}='IAAPloc2';
% exp(4)=2;
% 
% dimFileName='Discriminant_2spm_discriminantJL_spm_wavep_allm_06062019_EAPA.mat';
% try
%     dimPathName='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 5 Decision time-LDA all mouse/EAPA JL AM/';
%     load([dimPathName dimFileName])
% catch
%     dimPathName='C:\Users\Justin Losacco\Dropbox\Restrepo_Lab\Papers\JL paper\Identity paper\SPM manuscript  Justin-Daniel\Fig 5 Decision time-LDA all mouse\EAPA JL AM\';
%     load([dimPathName dimFileName])
% end
% odor_pair_no=odor_pair_no+1;
% handles_dim_per_odor_pair(odor_pair_no).handles_out=handles_out;
% odor_pair_label{5}='EAPAloc2';
% exp(5)=2;
% 
% dimFileName='Discriminant_2spm_discriminant_LFP_per_mouse_06102019_IsoAA_mo.mat';
% try
%     dimPathName='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 5 Decision time-LDA all mouse/IAMO JL AM/';
%     load([dimPathName dimFileName])
% catch
%     dimPathName='C:\Users\Justin Losacco\Dropbox\Restrepo_Lab\Papers\JL paper\Identity paper\SPM manuscript  Justin-Daniel\Fig 5 Decision time-LDA all mouse\IAMO JL AM\';
%     load([dimPathName dimFileName])
% end
% odor_pair_no=odor_pair_no+1;
% handles_dim_per_odor_pair(odor_pair_no).handles_out=handles_out;
% odor_pair_label{6}='IAMOloc2';
% exp(6)=2;

% PACnames{1}='Beta';
% PACnames{2}='Low gamma';
% PACnames{3}='High gamma';
% 
% window_label{1}='Pre-odor';
% window_label{2}='Odor';

prof_label{1}='Proficient';
prof_label{2}='Naive';


% t_odor_arrival=0.1;


%Used to troubleshoot the IAMO JL
% mice_included=[1 3 4];

figNo=0;
figNo1=12;
figNo2=24;
%Plot average percent correct for the LDA for peak and trough for
%wavelet power referenced to PAC phase
t=handles_out.t_power;
groupNo=1;
fprintf(1, ['\n\n'])

for PACii=[1 3]
    
    
    glm_ii=0;
    glm_dim=[];
    
    glm_peak_ii=0;
    glm_peak_dim=[];
  
    glm_trough_ii=0;
    glm_trough_dim=[];

    
    for ii_op=1:odor_pair_no
        
        ii_normpre_stats=0;
        p_normpre_dim_stats=[];
        
        %Plot decision time reationship
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        hFig=figure(figNo);
        set(hFig, 'units','normalized','position',[.4 .4 .25 .25])
        
        
        figNo1=figNo1+1;
        try
            close(figNo1)
        catch
        end
        hFig1=figure(figNo1);
        set(hFig1, 'units','normalized','position',[.4 .4 .25 .25])
        
        for percent_correct_ii=2:-1:1
            
           
            %Find the number of mice whose dimensionality was calculated for this odor
            %pair
            no_mice=0;
            mice_with_dim=[];
            dim_trough=[];
            dim_peak=[]; 
            for mouseNo=1:length(handles_dim_per_odor_pair(ii_op).handles_out.discriminant_PACwavepower)
                mice_with_dim(mouseNo)=0;
                if ~isempty(handles_dim_per_odor_pair(ii_op).handles_out.discriminant_PACpower_per_mouse(mouseNo).group)
                    try
                        if handles_dim_per_odor_pair(ii_op).handles_out.discriminant_PACpower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated==1
                            no_trials=handles_dim_per_odor_pair(ii_op).handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).no_trials;
                            if no_trials>=30
                                no_mice=no_mice+1;
                                mice_with_dim(mouseNo)=1;
                                dim_trough(no_mice,:)=handles_dim_per_odor_pair(ii_op).handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).dimensionality_trough;
                                dim_peak(no_mice,:)=handles_dim_per_odor_pair(ii_op).handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).dimensionality_peak;
                                
                            end
                        end
                    catch
                    end
                end
            end
            
            
            %Normalize the dimensionality
            for mouseNo=1:no_mice
                dim_trough(mouseNo,:)=dim_trough(mouseNo,:)/mean(dim_trough(mouseNo,t<0));
                dim_peak(mouseNo,:)=dim_peak(mouseNo,:)/mean(dim_peak(mouseNo,t<0));
            end
            
            %Plot using single mouse plots
            figure(figNo)
            subplot(1,2,3-percent_correct_ii)
            hold on
            for ii=1:size(dim_trough,1)
                plot(t,dim_trough(ii,:),'Color',[0.7 0.7 1])
            end
            plot(t,mean(dim_trough,1),'-b','LineWidth',2)
            
            for ii=1:size(dim_peak,1)
                plot(t,dim_peak(ii,:),'Color',[1 0.7 0.7])
            end
            plot(t,mean(dim_peak,1),'-r','LineWidth',2)
            
            plot([0 0],[0 2],'-k')
            odorhl=plot([0 2.5],[0.5 0.5],'-k','LineWidth',5);
            plot([2.5 2.5],[0 2],'-k')
            
            %             plot([t(1) t(end)],[no_mice no_mice],'-g')
            
            ylabel('Time (sec)')
            if percent_correct_ii==1
                legend('Trough','Peak')
            end
            ylim([0.5 2])
            
            
            
            %Plot using bounded lines
            figure(figNo1)
            subplot(1,2,3-percent_correct_ii)
            hold on
            %             for ii=1:size(dim_trough,1)
            %                 plot(t,dim_trough(ii,:),'Color',[0.7 0.7 1])
            %             end
            
            CItroughdims = bootci(1000, {@mean, dim_trough})';
            CItroughdims(:,1)=mean(dim_trough,1)'-CItroughdims(:,1);
            CItroughdims(:,2)=CItroughdims(:,2)-mean(dim_trough,1)';
            [hlCR, hpCR] = boundedline(t,mean(dim_trough,1)', CItroughdims, 'b');
            plot(t,mean(dim_trough,1),'-b','LineWidth',2)
            
            
            
            CIpeakdims = bootci(1000, {@mean, dim_peak})';
            CIpeakdims(:,1)=mean(dim_peak,1)'-CIpeakdims(:,1);
            CIpeakdims(:,2)=CIpeakdims(:,2)-mean(dim_peak,1)';
            [hlCR, hpCR] = boundedline(t,mean(dim_peak,1)', CIpeakdims, 'r');
            plot(t,mean(dim_peak,1),'-r','LineWidth',2)
            
            plot([0 0],[0 2],'-k')
            odorhl=plot([0 2.5],[0.5 0.5],'-k','LineWidth',5);
            plot([2.5 2.5],[0 2],'-k')
            
            %             plot([t(1) t(end)],[no_mice no_mice],'-g')
            
            ylabel('Time (sec)')
            if percent_correct_ii==1
                legend('Trough','Peak')
            end
            ylim([0.5 1.5])
            
            %Enter data for dimension glm
            dim_trough_pre=[];
            dim_peak_pre=[];
            dim_trough_odor=[];
            dim_peak_odor=[];
            
            for mouseNo=1:no_mice
                for ii_t=1:length(t)
                    %Trough
                    glm_ii=glm_ii+1;
                    glm_dim.data(glm_ii)=dim_trough(mouseNo,ii_t);
                    glm_dim.proficiency(glm_ii)=percent_correct_ii;
                    glm_dim.peak_trough(glm_ii)=0;
                    glm_dim.time(glm_ii)=t(ii_t);
                    
                    if ii_op<=3
                        glm_dim.exp(glm_ii)=1;
                    else
                        glm_dim.exp(glm_ii)=2;
                    end
                    
                    %Trough only
                    glm_trough_ii=glm_trough_ii+1;
                    glm_trough_dim.data(glm_trough_ii)=dim_trough(mouseNo,ii_t);
                    glm_trough_dim.proficiency(glm_trough_ii)=percent_correct_ii;
                    glm_trough_dim.time(glm_trough_ii)=t(ii_t);
                    
                    if ii_op<=3
                        glm_trough_dim.exp(glm_trough_ii)=1;
                    else
                        glm_trough_dim.exp(glm_trough_ii)=2;
                    end
                    
                    %Peak
                    glm_ii=glm_ii+1;
                    glm_dim.data(glm_ii)=dim_peak(mouseNo,ii_t);
                    glm_dim.proficiency(glm_ii)=percent_correct_ii;
                    glm_dim.peak_trough(glm_ii)=1;
                    glm_dim.time(glm_ii)=t(ii_t);
                    
                    if ii_op<=3
                        glm_dim.exp(glm_ii)=1;
                    else
                        glm_dim.exp(glm_ii)=2;
                    end
                    
                    %Peak only
                    glm_peak_ii=glm_peak_ii+1;
                    glm_peak_dim.data(glm_peak_ii)=dim_peak(mouseNo,ii_t);
                    glm_peak_dim.proficiency(glm_peak_ii)=percent_correct_ii;
                    glm_peak_dim.time(glm_peak_ii)=t(ii_t);
                    
                    if ii_op<=3
                        glm_peak_dim.exp(glm_peak_ii)=1;
                    else
                        glm_peak_dim.exp(glm_peak_ii)=2;
                    end
                end
                
                %I do the post hoc for two time periods: pre and odor
                dim_trough_pre(mouseNo)=mean(dim_trough(mouseNo,(t>=-1)&(t<0)),2);
                dim_peak_pre(mouseNo)=mean(dim_peak(mouseNo,(t>=-1)&(t<0)),2);
                dim_trough_odor(mouseNo)=mean(dim_trough(mouseNo,(t>=2)&(t<3)),2);
                dim_peak_odor(mouseNo)=mean(dim_peak(mouseNo,(t>=2)&(t<3)),2);  
            end
            
            %Save data for drgMutiRanksumorTtest for dimension normalized to pre-odor
            %I do the post hoc for two time points pre and odor
            
            %Pre-odor trough
            ii_normpre_stats=ii_normpre_stats+1;
            p_normpre_dim_stats(ii_normpre_stats).data=dim_trough_pre;
            p_normpre_dim_stats(ii_normpre_stats).description=['Trough pre-odor for ' odor_pair_label{ii_op} ' ' prof_label{percent_correct_ii}];
            
            %Pre-odor peak
            ii_normpre_stats=ii_normpre_stats+1;
            p_normpre_dim_stats(ii_normpre_stats).data=dim_peak_pre;
            p_normpre_dim_stats(ii_normpre_stats).description=['Peak pre-odor for ' odor_pair_label{ii_op} ' ' prof_label{percent_correct_ii}];
            
            %Odor trough
            ii_normpre_stats=ii_normpre_stats+1;
            p_normpre_dim_stats(ii_normpre_stats).data=dim_trough_odor;
            p_normpre_dim_stats(ii_normpre_stats).description=['Trough odor for ' odor_pair_label{ii_op} ' ' prof_label{percent_correct_ii}];
            
            %Odor peak
            ii_normpre_stats=ii_normpre_stats+1;
            p_normpre_dim_stats(ii_normpre_stats).data=dim_peak_odor;
            p_normpre_dim_stats(ii_normpre_stats).description=['Peak odor for ' odor_pair_label{ii_op} ' ' prof_label{percent_correct_ii}];
            
           fprintf(1, ['Number of mice for ' odor_pair_label{ii_op}  ' Theta/' handles_out.drgbchoices.PACnames{PACii} ' ' prof_label{percent_correct_ii} '= %d\n'], no_mice)
        end
        
        suptitle(['Norm dim for ' odor_pair_label{ii_op} 'Theta/' handles_out.drgbchoices.PACnames{PACii} ' per mouse'])
        
        %    Do post hoc tests
        fprintf(1, ['\n\nRanksum or t-test p values for dimensionality normalized to peak pre-odor for Theta/' PACnames{PACii} '\n'])
        try
            [output_data] = drgMutiRanksumorTtest(p_normpre_dim_stats);
            fprintf(1, '\n\n')
        catch
        end
    end
    
   
     
    %Perform the glm for dimensionality
    fprintf(1, ['\n\nglm for normalized dimensionality for Theta/' handles_out.drgbchoices.PACnames{PACii} '\n'])
    tbl = table(glm_dim.data',glm_dim.peak_trough',glm_dim.proficiency',glm_dim.exp',glm_dim.time',...
        'VariableNames',{'dimensionality','peak_trough','proficient_naive','exp','time'});
    mdl = fitglm(tbl,'dimensionality~peak_trough+proficient_naive+exp+time+peak_trough*proficient_naive*exp*time'...
        ,'CategoricalVars',[2,3,4])
    
    
   %Perform the glm for peak PRP dimensionality
    fprintf(1, ['\n\nglm for normalized peak PRP dimensionality for Theta/' handles_out.drgbchoices.PACnames{PACii} '\n'])
    tbl = table(glm_peak_dim.data',glm_peak_dim.proficiency',glm_peak_dim.exp',glm_peak_dim.time',...
        'VariableNames',{'dimensionality','proficient_naive','exp','time'});
    mdl = fitglm(tbl,'dimensionality~proficient_naive+exp+time+proficient_naive*exp*time'...
        ,'CategoricalVars',[2 3])

    
       %Perform the glm for trough PRP dimensionality
    fprintf(1, ['\n\nglm for normalized trough PRP dimensionality for Theta/' handles_out.drgbchoices.PACnames{PACii} '\n'])
    tbl = table(glm_trough_dim.data',glm_trough_dim.proficiency',glm_trough_dim.exp',glm_trough_dim.time',...
        'VariableNames',{'dimensionality','proficient_naive','exp','time'});
    mdl = fitglm(tbl,'dimensionality~proficient_naive+exp+time+proficient_naive*exp*time'...
        ,'CategoricalVars',[2,3])
    
    
end

pfff=1
   
    





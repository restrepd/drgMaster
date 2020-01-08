function new_all_mouse_dim_summary


warning('off')

close all
clear all

no_odor_pairs=6;

odor_pair_no=0;

%Daniel's
dimFileName='Discriminant_2spm_discriminantolfac_all_mice_aceto060519.mat';
try
    dimPathName='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 5 Decision time-LDA all mouse/APEB DR AM/';
    load([dimPathName dimFileName])
catch
    dimPathName='C:\Users\Justin Losacco\Dropbox\Restrepo_Lab\Papers\JL paper\Identity paper\SPM manuscript  Justin-Daniel\Fig 5 Decision time-LDA all mouse\APEB DR AM\';
    load([dimPathName dimFileName])
end
odor_pair_no=odor_pair_no+1;
handles_dim_per_odor_pair(odor_pair_no).handles_out=handles_out;
odor_pair_label{1}='APEBloc1';
exp(1)=1;

dimFileName='Discriminant_2spmc_discriminantolfac_all_mice_EAPA06062019.mat';
try
    dimPathName='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 5 Decision time-LDA all mouse/EAPA DR AM/';
    load([dimPathName dimFileName])
catch
    dimPathName='C:\Users\Justin Losacco\Dropbox\Restrepo_Lab\Papers\JL paper\Identity paper\SPM manuscript  Justin-Daniel\Fig 5 Decision time-LDA all mouse\EAPA DR AM\';
    load([dimPathName dimFileName])
end
odor_pair_no=odor_pair_no+1;
handles_dim_per_odor_pair(odor_pair_no).handles_out=handles_out;
odor_pair_label{2}='EAPAloc1';
exp(2)=1;

dimFileName='Discriminant_2spm_discriminantolfac_all_mice_iso06082019.mat';
try
    dimPathName='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 5 Decision time-LDA all mouse/IAMO DR AM/';
    load([dimPathName dimFileName])
catch
    dimPathName='C:\Users\Justin Losacco\Dropbox\Restrepo_Lab\Papers\JL paper\Identity paper\SPM manuscript  Justin-Daniel\Fig 5 Decision time-LDA all mouse\IAMO DR AM\';
    load([dimPathName dimFileName])
end
odor_pair_no=odor_pair_no+1;
handles_dim_per_odor_pair(odor_pair_no).handles_out=handles_out;
odor_pair_label{3}='IAMOloc1';
exp(3)=1;

%Justin's
dimFileName='Discriminant_2spm_discriminant_LFP_wp_all_mouse_05312019_IAAP.mat';
try
    dimPathName='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 5 Decision time-LDA all mouse/IAAP JL AM/';
    load([dimPathName dimFileName])
catch
    dimPathName='C:\Users\Justin Losacco\Dropbox\Restrepo_Lab\Papers\JL paper\Identity paper\SPM manuscript  Justin-Daniel\Fig 5 Decision time-LDA all mouse\IAAP JL AM\';
    load([dimPathName dimFileName])
end
odor_pair_no=odor_pair_no+1;
handles_dim_per_odor_pair(odor_pair_no).handles_out=handles_out;
odor_pair_label{4}='IAAPloc2';
exp(4)=2;

dimFileName='Discriminant_2spm_discriminantJL_spm_wavep_allm_06062019_EAPA.mat';
try
    dimPathName='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 5 Decision time-LDA all mouse/EAPA JL AM/';
    load([dimPathName dimFileName])
catch
    dimPathName='C:\Users\Justin Losacco\Dropbox\Restrepo_Lab\Papers\JL paper\Identity paper\SPM manuscript  Justin-Daniel\Fig 5 Decision time-LDA all mouse\EAPA JL AM\';
    load([dimPathName dimFileName])
end
odor_pair_no=odor_pair_no+1;
handles_dim_per_odor_pair(odor_pair_no).handles_out=handles_out;
odor_pair_label{5}='EAPAloc2';
exp(5)=2;

dimFileName='Discriminant_2spm_discriminant_LFP_all_mouse_06102019_IsoAA_mo.mat';
try
    dimPathName='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 5 Decision time-LDA all mouse/IAMO JL AM/';
    load([dimPathName dimFileName])
catch
    dimPathName='C:\Users\Justin Losacco\Dropbox\Restrepo_Lab\Papers\JL paper\Identity paper\SPM manuscript  Justin-Daniel\Fig 5 Decision time-LDA all mouse\IAMO JL AM\';
    load([dimPathName dimFileName])
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

prof_label{1}='Proficient';
prof_label{2}='Naive'


t_odor_arrival=0.1;


figNo1=0;
figNo2=4;

%Plot average percent correct for the LDA for peak and trough for
%wavelet power referenced to PAC phase
t=handles_out.t_power;
groupNo=1;

for PACii=[1 3]
    
    
    
    all_dim_trough=zeros(2,2,3,71);
    all_dim_peak=zeros(2,2,3,71);
    all_normmousedim_trough=zeros(2,2,3,71);
    all_normmousedim_peak=zeros(2,2,3,71);
    all_normpredim_trough=zeros(2,2,3,71);
    all_normpredim_peak=zeros(2,2,3,71);
    ii_loc1=zeros(1,2);
    ii_loc2=zeros(1,2);
    
%     glm_ii=0;
%     glm_dim=[];
%     ii_dim_stats=0;
%     p_dim_stats=[];
%     
%     glm_norm_ii=0;
%     glm_norm_dim=[];
%     ii_norm_stats=0;
%     p_norm_dim_stats=[];
    
    glm_normpre_ii=0;
    glm_normpre_dim=[];
    
    glm_peak_ii=0;
    glm_peak_dim=[];
    
    glm_trough_ii=0;
    glm_trough_dim=[];
    
%     ii_normpre_stats=0;
%     p_normpre_dim_stats=[];
%     p_normpre_dim_stats_all_exp=[];
    
    for ii_op=1:odor_pair_no
        
        
        %         %Plot decision time reationship
        %         figNo=figNo+1;
        %         try
        %             close(figNo)
        %         catch
        %         end
        %         hFig=figure(figNo);
        %         set(hFig, 'units','normalized','position',[.4 .4 .25 .25])
        %         hold on
        
        for percent_correct_ii=2:-1:1
            
            
            %Find the number of mice
            no_mice=0;
            for mouseNo=1:length(handles_out.all_mouse_wav)
                try
                    no_trials=handles_out.all_mouse_wav(mouseNo).group(groupNo).percent_correct(percent_correct_ii).no_trials;
                    if no_trials>=30
                        no_mice=no_mice+1;
                    end
                catch
                end
            end
            
            dim_trough=handles_dim_per_odor_pair(ii_op).handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).dimensionality_trough;
            dim_peak=handles_dim_per_odor_pair(ii_op).handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).dimensionality_peak;
            dim_peak_proficient=handles_dim_per_odor_pair(ii_op).handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(1).PACii(PACii).dimensionality_peak;
            no_mice=handles_dim_per_odor_pair(ii_op).handles_out.discriminant_PACwavepower_all_mice.no_mice_included;
            
            if ii_op<=3
                ii_loc1(percent_correct_ii)=ii_loc1(percent_correct_ii)+1;
                %                 all_dim_trough(percent_correct_ii,1,ii_loc1(percent_correct_ii),:)=dim_trough;
                %                 all_dim_peak(percent_correct_ii,1,ii_loc1(percent_correct_ii),:)=dim_peak;
                %                 all_normmousedim_trough(percent_correct_ii,1,ii_loc1(percent_correct_ii),:)=dim_trough/no_mice;
                %                 all_normmousedim_peak(percent_correct_ii,1,ii_loc1(percent_correct_ii),:)=dim_peak/no_mice;
                %                 all_normpredim_trough(percent_correct_ii,1,ii_loc1(percent_correct_ii),:)=dim_trough/mean(dim_peak_proficient(t<0));
                all_normpredim_trough(percent_correct_ii,1,ii_loc1(percent_correct_ii),:)=dim_trough/mean(dim_trough(t<0));
                %                 all_normpredim_peak(percent_correct_ii,1,ii_loc1(percent_correct_ii),:)=dim_peak/mean(dim_peak_proficient(t<0));
                all_normpredim_peak(percent_correct_ii,1,ii_loc1(percent_correct_ii),:)=dim_peak/mean(dim_peak(t<0));
            else
                ii_loc2(percent_correct_ii)=ii_loc2(percent_correct_ii)+1;
                %                 all_dim_trough(percent_correct_ii,2,ii_loc2(percent_correct_ii),:)=dim_trough;
                %                 all_dim_peak(percent_correct_ii,2,ii_loc2(percent_correct_ii),:)=dim_peak;
                %                 all_normmousedim_trough(percent_correct_ii,2,ii_loc2(percent_correct_ii),:)=dim_trough/no_mice;
                %                 all_normmousedim_peak(percent_correct_ii,2,ii_loc2(percent_correct_ii),:)=dim_peak/no_mice;
                %                 all_normpredim_trough(percent_correct_ii,2,ii_loc2(percent_correct_ii),:)=dim_trough/mean(dim_peak_proficient(t<0));
                %                 all_normpredim_peak(percent_correct_ii,2,ii_loc2(percent_correct_ii),:)=dim_peak/mean(dim_peak_proficient(t<0));
                all_normpredim_trough(percent_correct_ii,2,ii_loc2(percent_correct_ii),:)=dim_trough/mean(dim_trough(t<0));
                all_normpredim_peak(percent_correct_ii,2,ii_loc2(percent_correct_ii),:)=dim_peak/mean(dim_peak(t<0));
            end
            
            %             subplot(1,2,3-percent_correct_ii)
            %             hold on
            %             plot(t,dim_trough,'-b')
            %             plot(t,dim_peak,'-r')
            %
            %             plot([0 0],[0 18],'-k')
            %             odorhl=plot([0 2.5],[0.5 0.5],'-k','LineWidth',5);
            %             plot([2.5 2.5],[0 18],'-k')
            %
            % %             plot([t(1) t(end)],[no_mice no_mice],'-g')
            %
            %             ylabel('Time (sec)')
            %             if percent_correct_ii==1
            %                 legend('Trough','Peak')
            %             end
            %             ylim([0 18])
            
            %             %Enter data for dimension glm
            %             for ii_t=1:length(t)
            %                 %Trough
            %                 glm_ii=glm_ii+1;
            %                 glm_dim.data(glm_ii)=dim_trough(ii_t);
            %                 glm_dim.proficiency(glm_ii)=percent_correct_ii;
            %                 glm_dim.peak_trough(glm_ii)=0;
            %                 glm_dim.time(glm_ii)=t(ii_t);
            %
            %                 if ii_op<=3
            %                     glm_dim.exp(glm_ii)=1;
            %                 else
            %                     glm_dim.exp(glm_ii)=2;
            %                 end
            %
            %                 %Peak
            %                 glm_ii=glm_ii+1;
            %                 glm_dim.data(glm_ii)=dim_peak(ii_t);
            %                 glm_dim.proficiency(glm_ii)=percent_correct_ii;
            %                 glm_dim.peak_trough(glm_ii)=1;
            %                 glm_dim.time(glm_ii)=t(ii_t);
            %
            %                 if ii_op<=3
            %                     glm_dim.exp(glm_ii)=1;
            %                 else
            %                     glm_dim.exp(glm_ii)=2;
            %                 end
            %             end
            
            %Enter data for glm for dimension normalized to pre-odor
            for ii_t=1:length(t)
                %Trough
                glm_normpre_ii=glm_normpre_ii+1;
                glm_normpre_dim.data(glm_normpre_ii)=dim_trough(ii_t)/mean(dim_trough(t<0));
                glm_normpre_dim.proficiency(glm_normpre_ii)=percent_correct_ii;
                glm_normpre_dim.peak_trough(glm_normpre_ii)=0;
                glm_normpre_dim.time(glm_normpre_ii)=t(ii_t);
                
                if ii_op<=3
                    glm_normpre_dim.exp(glm_normpre_ii)=1;
                else
                    glm_normpre_dim.exp(glm_normpre_ii)=2;
                end
                
                %Only trough
                glm_trough_ii=glm_trough_ii+1;
                glm_trough_dim.data(glm_trough_ii)=dim_trough(ii_t)/mean(dim_trough(t<0));
                glm_trough_dim.proficiency(glm_trough_ii)=percent_correct_ii;
                glm_trough_dim.time(glm_trough_ii)=t(ii_t);
                
                if ii_op<=3
                    glm_trough_dim.exp(glm_trough_ii)=1;
                else
                    glm_trough_dim.exp(glm_trough_ii)=2;
                end
                
                %Peak
                glm_normpre_ii=glm_normpre_ii+1;
                glm_normpre_dim.data(glm_normpre_ii)=dim_peak(ii_t)/mean(dim_peak(t<0));
                glm_normpre_dim.proficiency(glm_normpre_ii)=percent_correct_ii;
                glm_normpre_dim.peak_trough(glm_normpre_ii)=1;
                glm_normpre_dim.time(glm_normpre_ii)=t(ii_t);
                
                if ii_op<=3
                    glm_normpre_dim.exp(glm_normpre_ii)=1;
                else
                    glm_normpre_dim.exp(glm_normpre_ii)=2;
                end
                
                %Only peak
                glm_peak_ii=glm_peak_ii+1;
                glm_peak_dim.data(glm_peak_ii)=dim_peak(ii_t)/mean(dim_peak(t<0));
                glm_peak_dim.proficiency(glm_peak_ii)=percent_correct_ii;
                glm_peak_dim.time(glm_peak_ii)=t(ii_t);
                
                if ii_op<=3
                    glm_peak_dim.exp(glm_peak_ii)=1;
                else
                    glm_peak_dim.exp(glm_peak_ii)=2;
                end
                
            end
            
            
            %Save data for drgMutiRanksumorTtest for dimension normalized to pre-odor
            
            %Pre-odor trough
            dim_trough_prenorm=dim_trough/mean(dim_trough(t<0));
            p_normpre_dim_stats(4*(percent_correct_ii-1)+1).data(ii_op)=mean(dim_trough_prenorm((t>=-1)&(t<0)));
            p_normpre_dim_stats(4*(percent_correct_ii-1)+1).description=['Trough pre-odor ' prof_label{percent_correct_ii}];
            
            
            %Pre-odor peak
            dim_peak_prenorm=dim_peak/mean(dim_peak(t<0));
            p_normpre_dim_stats(4*(percent_correct_ii-1)+2).data(ii_op)=mean(dim_peak_prenorm((t>=-1)&(t<0)));
            p_normpre_dim_stats(4*(percent_correct_ii-1)+2).description=['Peak pre-odor ' prof_label{percent_correct_ii}];
            
            
            %Odor trough
            p_normpre_dim_stats(4*(percent_correct_ii-1)+3).data(ii_op)=mean(dim_trough_prenorm((t>=2)&(t<3)));
            p_normpre_dim_stats(4*(percent_correct_ii-1)+3).description=['Trough odor ' prof_label{percent_correct_ii}];
            
            
            %Odor peak
            p_normpre_dim_stats(4*(percent_correct_ii-1)+4).data(ii_op)=mean(dim_peak_prenorm((t>=2)&(t<3)));
            p_normpre_dim_stats(4*(percent_correct_ii-1)+4).description=['Peak odor ' prof_label{percent_correct_ii}];
            
            
        end
        
        %     suptitle(['Dimensionality for ' odor_pair_label{ii_op} 'Theta/' handles_out.drgbchoices.PACnames{PACii} ' all mice pooled'])
    end
    
    
    %Perform the glm for dimensionality normalized to peal pre-odor
    fprintf(1, ['\n\nglm for normalized dimensionality for Theta/' handles_out.drgbchoices.PACnames{PACii} '\n'])
    tbl = table(glm_normpre_dim.data',glm_normpre_dim.peak_trough',glm_normpre_dim.proficiency',glm_normpre_dim.exp',glm_normpre_dim.time',...
        'VariableNames',{'dimensionality','peak_trough','proficient_naive','exp','time'});
    mdl = fitglm(tbl,'dimensionality~peak_trough+proficient_naive+exp+time+peak_trough*proficient_naive*exp*time'...
        ,'CategoricalVars',[2,3,4])
    
      %Perform the glm for peak normalized dimensionality
    fprintf(1, ['\n\nglm for normalized dimensionality for peak for Theta/' handles_out.drgbchoices.PACnames{PACii} '\n'])
    tbl = table(glm_peak_dim.data',glm_peak_dim.proficiency',glm_peak_dim.exp',glm_peak_dim.time',...
        'VariableNames',{'dimensionality','proficient_naive','exp','time'});
    mdl = fitglm(tbl,'dimensionality~proficient_naive+exp+time+proficient_naive*exp*time'...
        ,'CategoricalVars',[2,3])
    
    %Perform the glm for trough normalized dimensionality
    fprintf(1, ['\n\nglm for normalized dimensionality for trough for Theta/' handles_out.drgbchoices.PACnames{PACii} '\n'])
    tbl = table(glm_trough_dim.data',glm_trough_dim.proficiency',glm_trough_dim.exp',glm_trough_dim.time',...
        'VariableNames',{'dimensionality','proficient_naive','exp','time'});
    mdl = fitglm(tbl,'dimensionality~proficient_naive+exp+time+proficient_naive*exp*time'...
        ,'CategoricalVars',[2,3])
    
    %Note: No post-hoc t-test/ranksum becuse there are only three points for
    %each experimental set
    %         Do ranksum/t test
    fprintf(1, ['\n\nRanksum or t-test p values for dimensionality normalized to peak pre-odor for both experiments for Theta/' PACnames{PACii} '\n'])
    try
        [output_data] = drgMutiRanksumorTtest(p_normpre_dim_stats);
        fprintf(1, '\n\n')
    catch
    end
    
    
    pffft=1;
    
    
    
    %Plot decision time reationship for normalized dimension averages for exp1
    %and exp2 including plots for each odorant
    
    for exp=1:2
        figNo1=figNo1+1;
        try
            close(figNo1)
        catch
        end
        hFig=figure(figNo1);
        set(hFig, 'units','normalized','position',[.4 .4 .25 .25])
        hold on
        
        for percent_correct_ii=2:-1:1
            
            subplot(1,2,3-percent_correct_ii)
            hold on
            
            
            these_trough_dims=zeros(3,71);
            these_trough_dims(:,:)=all_normpredim_trough(percent_correct_ii,exp,:,:);
            mean_trouhg_dims=mean(these_trough_dims)';
            %             CItroughdims = bootci(1000, {@mean, these_trough_dims})';
            %             CItroughdims(:,1)=mean_trouhg_dims-CItroughdims(:,1);
            %             CItroughdims(:,2)=CItroughdims(:,2)-mean_trouhg_dims;
            %             [hlCR, hpCR] = boundedline(t,mean_trouhg_dims, CItroughdims, 'b');
            p1=plot(t,mean_trouhg_dims, 'b','LineWidth',2);
            for ii_exp=1:size(these_trough_dims,1)
                plot(t,these_trough_dims(ii_exp,:),'Color',[0.7 0.7 1])
            end
            
            
            these_peak_dims=zeros(3,71);
            these_peak_dims(:,:)=all_normpredim_peak(percent_correct_ii,exp,:,:);
            mean_peak_dims=mean(these_peak_dims)';
            %             CIpeakdims = bootci(1000, {@mean, these_peak_dims})';
            %             CIpeakdims(:,1)=mean_peak_dims-CIpeakdims(:,1);
            %             CIpeakdims(:,2)=CIpeakdims(:,2)-mean_peak_dims;
            %             [hlCR, hpCR] = boundedline(t,mean_peak_dims, CIpeakdims, 'r');
            p2=plot(t,mean_peak_dims, 'r','LineWidth',2);
            for ii_exp=1:size(these_peak_dims,1)
                plot(t,these_peak_dims(ii_exp,:),'Color',[1 0.7 0.7])
            end
            
            
            plot([0 0],[0 18],'-k')
            odorhl=plot([0 2.5],[0.2 0.2],'-k','LineWidth',5);
            plot([2.5 2.5],[0 18],'-k')
            
            %         plot([t(1) t(end)],[no_mice no_mice],'-g')
            
            ylabel('Time (sec)')
            if percent_correct_ii==1
                legend([p1 p2], {'Trough','Peak'})
            end
            ylim([0 1.5])
            
        end
        if exp==1
            suptitle(['Dim-pooled-mice/pre-dim for Theta/' handles_out.drgbchoices.PACnames{PACii} ' exp1'])
        else
            suptitle(['Dim-pooled-mince/pre-dim for Theta/' handles_out.drgbchoices.PACnames{PACii} ' exp2'])
        end
    end
    
    
    %Plot decision time reationship for normalized dimension averages for exp1
    %and exp2 using bounded lines
    
    for exp=1:2
        figNo2=figNo2+1;
        try
            close(figNo2)
        catch
        end
        hFig=figure(figNo2);
        set(hFig, 'units','normalized','position',[.4 .4 .25 .25])
        hold on
        
        for percent_correct_ii=2:-1:1
            
            subplot(1,2,3-percent_correct_ii)
            hold on
            
            
            these_trough_dims=zeros(3,71);
            these_trough_dims(:,:)=all_normpredim_trough(percent_correct_ii,exp,:,:);
            mean_trouhg_dims=mean(these_trough_dims)';
            CItroughdims = bootci(1000, {@mean, these_trough_dims})';
            CItroughdims(:,1)=mean_trouhg_dims-CItroughdims(:,1);
            CItroughdims(:,2)=CItroughdims(:,2)-mean_trouhg_dims;
            [hlCR, hpCR] = boundedline(t,mean_trouhg_dims, CItroughdims, 'b');
            p1=plot(t,mean_trouhg_dims, 'b','LineWidth',2);
            %             for ii_exp=1:size(these_trough_dims,1)
            %                 plot(t,these_trough_dims(ii_exp,:),'Color',[0.7 0.7 1])
            %             end
            %
            
            these_peak_dims=zeros(3,71);
            these_peak_dims(:,:)=all_normpredim_peak(percent_correct_ii,exp,:,:);
            mean_peak_dims=mean(these_peak_dims)';
            CIpeakdims = bootci(1000, {@mean, these_peak_dims})';
            CIpeakdims(:,1)=mean_peak_dims-CIpeakdims(:,1);
            CIpeakdims(:,2)=CIpeakdims(:,2)-mean_peak_dims;
            [hlCR, hpCR] = boundedline(t,mean_peak_dims, CIpeakdims, 'r');
            p2=plot(t,mean_peak_dims, 'r','LineWidth',2);
            %              for ii_exp=1:size(these_peak_dims,1)
            %                 plot(t,these_peak_dims(ii_exp,:),'Color',[1 0.7 0.7])
            %             end
            
            
            plot([0 0],[0 18],'-k')
            odorhl=plot([0 2.5],[0.2 0.2],'-k','LineWidth',5);
            plot([2.5 2.5],[0 18],'-k')
            
            %         plot([t(1) t(end)],[no_mice no_mice],'-g')
            
            ylabel('Time (sec)')
            if percent_correct_ii==1
                legend([p1 p2], {'Trough','Peak'})
            end
            ylim([0 1.5])
            
        end
        if exp==1
            suptitle(['Dim-pooled-mice/pre-dim for Theta/' handles_out.drgbchoices.PACnames{PACii} ' exp1'])
        else
            suptitle(['Dim-pooled-mince/pre-dim for Theta/' handles_out.drgbchoices.PACnames{PACii} ' exp2'])
        end
    end
    
end



   
    





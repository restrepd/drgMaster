function drgAnalyzeLFPDiscriminantMultiBatch
%Analyzes the linear discriminant analysis performed by drgLFPDiscriminantBatch
%Takes as a in input the 'drgbDiscPar' file listing 'Discriminant_*.mat' output files
%from drgLFPDiscriminantBatch
%
%Performs summary analises for LDA and  PCA

% which_display chooses the analysis:
%
%1 Displays average predicton for proficeint vs naive for LDA and PCA for power LFP
%
%2 Displays average predicton for proficeint vs naive for LDA and PCA for angle in PAC
%
%3 Displays average prediction and dimensionality for peak and trough for LDA for wavelet
%power referenced to the phase of PAC and plots PC1 for the PCA.
%These are choices 10 and 11 in drgLFPDiscriminantBatch

close all
clear all


which_display=3;
mice_excluded=[];

[choiceFileName,choiceBatchPathName] = uigetfile({'drgbDiscPar*.m'},'Select the .m file with all the choices for analysis');
fprintf(1, ['\ndrgAnalyzeLFPDiscriminantMultiBatch run for ' choiceFileName '\n\n']);

tempDirName=['temp' choiceFileName(12:end-2)];

addpath(choiceBatchPathName)
eval(['handlesdrgb=' choiceFileName(1:end-2) ';'])

groupNo=1; %Note I am doing only forward here

%Define the windows for analysis of dimensionality
window_start=[-1 0.5];
window_end=[0 2.5];
no_wins=2;

window_legends{1}='Pre-odor';
window_legends{2}='Odor';

%This is the window for area under the curve case 3
auc_from=0.1;
auc_to=2.5;

these_colors{1}='b';
these_colors{2}='r';
these_colors{3}='m';
these_colors{4}='g';
these_colors{5}='y';
these_colors{6}='k';
these_colors{7}='c';

%Plot the summary wavelet power at peak vs power at trough
for figNo=1:6
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    set(hFig, 'units','normalized','position',[.3 .3 .3 .3])
    hold on
end




for fileNo=1:handlesdrgb.drgbchoices.no_files
    
    pname=handlesdrgb.drgbchoices.PathName{fileNo};
    fname=handlesdrgb.drgbchoices.FileName{fileNo};
    
    discriminant_name=[pname fname];
    load(discriminant_name)
    
    
    
    
    
    %Plot average percent correct for the LDA for peak and trough for
    %wavelet power referenced to PAC phase
    t=handles_out.t_power;
    
    for PACii=1:length(handles_out.drgbchoices.PACburstLowF)
        p_correct_stats=[];
        ii_stats=0;
        p_dim_stats=[];
        ii_dim_stats=0;
        glm_ii=0;
        glm_correct=[];
        glm_dim_ii=0;
        glm_dim=[];
        for percent_correct_ii=1:2
            
            figure((PACii-1)*2+percent_correct_ii)
            
            %Gather all the data
            no_mice=0;
            no_mice_included=0;
            all_discriminant_correct_peak=zeros(max(handles_out.drgbchoices.mouse_no),length(t));
            all_discriminant_correct_trough=zeros(max(handles_out.drgbchoices.mouse_no),length(t));
            all_discriminant_correct_shuffled_peak=zeros(max(handles_out.drgbchoices.mouse_no),length(t));
            all_discriminant_correct_shuffled_trough=zeros(max(handles_out.drgbchoices.mouse_no),length(t));
        
            
            for mouseNo=1:length(handles_out.discriminant_PACwavepower)
                try
                    if handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated==1
                        per_ii=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).no_trials;
                        no_mice=no_mice+1;
                        if (per_ii>=20)&(sum(no_mice==mice_excluded)==0)
                            no_mice_included=no_mice_included+1;
                            all_discriminant_correct_peak(no_mice_included,:)=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_peak;
                            all_discriminant_correct_trough(no_mice_included,:)=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_trough;
                            all_discriminant_correct_shuffled_peak(no_mice_included,:)=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_shuffled_peak;
                            all_discriminant_correct_shuffled_trough(no_mice_included,:)=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_shuffled_trough;
                            
                        end
                    end
                catch
                end
            end
            
            
            no_mice_per(percent_correct_ii)=no_mice_included;
            
            %Plot percent correct for the LDA and save the data for
            %the ranksum
            %                 figNo=figNo+1;
            %                 try
            %                     close(figNo)
            %                 catch
            %                 end
            %                 hFig=figure(figNo);
            %
            %                 hold on
            
            %Note that I merge the shuffled for both peak and
            %trough for the average plot
            all_discriminant_correct_shuffled=zeros(2*no_mice_included,length(t));
            all_discriminant_correct_shuffled(1:no_mice_included,:)=all_discriminant_correct_shuffled_peak(1:no_mice_included,:);
            all_discriminant_correct_shuffled(no_mice_included+1:end,:)=all_discriminant_correct_shuffled_trough(1:no_mice_included,:);
            mean_dcsh=mean(all_discriminant_correct_shuffled,1)';
            %                 if size(all_discriminant_correct_shuffled,1)>2
            %                     CIdcsh = bootci(1000, {@mean, all_discriminant_correct_shuffled})';
            %                     CIdcsh(:,1)=mean_dcsh-CIdcsh(:,1);
            %                     CIdcsh(:,2)=CIdcsh(:,2)-mean_dcsh;
            %                     [hlCR, hpCR] = boundedline(t,mean_dcsh, CIdcsh, 'k');
            %                 else
            %                     plot(t,mean_dcsh,'-k')
            %                 end
            
            data=[];
            for mouseNo=1:no_mice_included
                data=[data (mean(all_discriminant_correct_shuffled_peak(mouseNo,(t>=auc_from)&(t<=auc_to)),2)-50)/50];
            end
            
            ii_stats=ii_stats+1;
            p_correct_stats(ii_stats).data=data;
            p_correct_stats(ii_stats).description=[handles_out.drgbchoices.group_no_names{groupNo}  ' ' handles_out.drgbchoices.per_lab{percent_correct_ii} ' shuffled peak'];
            p_correct_stats(ii_stats).data_ii=1;
            p_correct_stats(ii_stats).PACii=PACii;
            p_correct_stats(ii_stats).per_corr_ii=percent_correct_ii;
            p_correct_stats(ii_stats).groupNo=groupNo;
            %
            %
            %                 glm_correct.data(glm_ii+1:glm_ii+length(data))=data;
            %                 glm_correct.PACii(glm_ii+1:glm_ii+length(data))=PACii;
            %                 glm_correct.group(glm_ii+1:glm_ii+length(data))=groupNo;
            %                 glm_correct.perCorr(glm_ii+1:glm_ii+length(data))=percent_correct_ii;
            %                 glm_correct.shuffled(glm_ii+1:glm_ii+length(data))=1;
            %                 glm_correct.peak(glm_ii+1:glm_ii+length(data))=1;
            %                 glm_ii=glm_ii+length(data);
            %
            data=[];
            for mouseNo=1:no_mice_included
                data=[data (mean(all_discriminant_correct_shuffled_trough(mouseNo,(t>=auc_from)&(t<=auc_to)),2)-50)/50];
            end
            
            ii_stats=ii_stats+1;
            p_correct_stats(ii_stats).data=data;
            p_correct_stats(ii_stats).description=[handles_out.drgbchoices.group_no_names{groupNo}  ' ' handles_out.drgbchoices.per_lab{percent_correct_ii} ' shuffled trough'];
            p_correct_stats(ii_stats).data_ii=2;
            p_correct_stats(ii_stats).PACii=PACii;
            p_correct_stats(ii_stats).per_corr_ii=percent_correct_ii;
            p_correct_stats(ii_stats).groupNo=groupNo;
            %
            %                 glm_correct.data(glm_ii+1:glm_ii+length(data))=data;
            %                 glm_correct.PACii(glm_ii+1:glm_ii+length(data))=PACii;
            %                 glm_correct.group(glm_ii+1:glm_ii+length(data))=groupNo;
            %                 glm_correct.perCorr(glm_ii+1:glm_ii+length(data))=percent_correct_ii;
            %                 glm_correct.shuffled(glm_ii+1:glm_ii+length(data))=1;
            %                 glm_correct.peak(glm_ii+1:glm_ii+length(data))=0;
            %                 glm_ii=glm_ii+length(data);
            %
            %
            %
            %Now plot the percent correct for the trough
            all_discriminant_correct_trough=all_discriminant_correct_trough(1:no_mice_included,:);
            mean_dc_trough=mean(all_discriminant_correct_trough,1)';
            %                 if size(all_discriminant_correct_trough,1)>2
            %                     CIdc_trough = bootci(1000, {@mean, all_discriminant_correct_trough})';
            %                     CIdc_trough(:,1)=mean_dc_trough-CIdc_trough(:,1);
            %                     CIdc_trough(:,2)=CIdc_trough(:,2)-mean_dc_trough;
            %                     [hlCR, hpCR] = boundedline(t,mean_dc_trough, CIdc_trough, 'b');
            %                 else
            %                     plot(t,mean_dc_trough,'r')
            %                 end
            %
            %
            data=[];
            for mouseNo=1:no_mice_included
                data=[data (mean(all_discriminant_correct_trough(mouseNo,(t>=auc_from)&(t<=auc_to)),2)-50)/50];
            end
            ii_stats=ii_stats+1;
            p_correct_stats(ii_stats).data=data;
            p_correct_stats(ii_stats).description=[handles_out.drgbchoices.group_no_names{groupNo} ' ' handles_out.drgbchoices.per_lab{percent_correct_ii} ' trough'];
            p_correct_stats(ii_stats).data_ii=3;
            p_correct_stats(ii_stats).PACii=PACii;
            p_correct_stats(ii_stats).per_corr_ii=percent_correct_ii;
            p_correct_stats(ii_stats).groupNo=groupNo;
            %
            %                 glm_correct.data(glm_ii+1:glm_ii+length(data))=data;
            %                 glm_correct.PACii(glm_ii+1:glm_ii+length(data))=PACii;
            %                 glm_correct.group(glm_ii+1:glm_ii+length(data))=groupNo;
            %                 glm_correct.perCorr(glm_ii+1:glm_ii+length(data))=percent_correct_ii;
            %                 glm_correct.shuffled(glm_ii+1:glm_ii+length(data))=0;
            %                 glm_correct.peak(glm_ii+1:glm_ii+length(data))=0;
            %                 glm_ii=glm_ii+length(data);
            
            %Now plot the percent correct for the peak
            all_discriminant_correct_peak=all_discriminant_correct_peak(1:no_mice_included,:);
            mean_dc_peak=mean(all_discriminant_correct_peak,1)';
            %                 if size(all_discriminant_correct_peak,1)>2
            %                     CIdc_peak = bootci(1000, {@mean, all_discriminant_correct_peak})';
            %                     CIdc_peak(:,1)=mean_dc_peak-CIdc_peak(:,1);
            %                     CIdc_peak(:,2)=CIdc_peak(:,2)-mean_dc_peak;
            %                     [hlCR, hpCR] = boundedline(t,mean_dc_peak, CIdc_peak, 'r');
            %                 else
            %                     plot(t,mean_dc_peak,'r')
            %                 end
            
            
            data=[];
            for mouseNo=1:no_mice_included
                data=[data (mean(all_discriminant_correct_peak(mouseNo,(t>=auc_from)&(t<=auc_to)),2)-50)/50];
            end
            ii_stats=ii_stats+1;
            p_correct_stats(ii_stats).data=data;
            p_correct_stats(ii_stats).description=[handles_out.drgbchoices.group_no_names{groupNo} ' ' ...
                handles_out.drgbchoices.per_lab{percent_correct_ii}...
                ' peak'];
            p_correct_stats(ii_stats).data_ii=4;
            p_correct_stats(ii_stats).PACii=PACii;
            p_correct_stats(ii_stats).per_corr_ii=percent_correct_ii;
            p_correct_stats(ii_stats).groupNo=groupNo;
            %
            %                 glm_correct.data(glm_ii+1:glm_ii+length(data))=data;
            %                 glm_correct.PACii(glm_ii+1:glm_ii+length(data))=PACii;
            %                 glm_correct.group(glm_ii+1:glm_ii+length(data))=groupNo;
            %                 glm_correct.perCorr(glm_ii+1:glm_ii+length(data))=percent_correct_ii;
            %                 glm_correct.shuffled(glm_ii+1:glm_ii+length(data))=0;
            %                 glm_correct.peak(glm_ii+1:glm_ii+length(data))=1;
            %                 glm_ii=glm_ii+length(data);
            
            %                 %Odor on markers
            %                 plot([0 0],[0 100],'-k')
            %                 odorhl=plot([0 2.5],[10 10],'-k','LineWidth',5);
            %                 plot([2.5 2.5],[0 100],'-k')
            %
            %                 title(['LDA for teta/' handles_out.drgbchoices.PACnames{PACii} ' ' handles_out.drgbchoices.group_no_names{groupNo}  ' ' handles_out.drgbchoices.per_lab{percent_correct_ii}])
            %
            %                 xlabel('Time (sec)')
            %                 ylabel(['% correct '  handles_out.drgbchoices.per_lab{percent_correct_ii}])
            %                 legend('Shuffled','Trough','Peak')
            
            %Plot the peak vs trough graph
            %                 figNo=figNo+1;
            %                 try
            %                     close(figNo)
            %                 catch
            %                 end
            %                 hFig=figure(figNo);
            
            hold on
            
            plot(p_correct_stats(ii_stats-1).data,p_correct_stats(ii_stats).data,'o','MarkerEdgeColor',these_colors{fileNo},'MarkerFaceColor',these_colors{fileNo})
            %plot(p_correct_stats(ii_stats-2).data,p_correct_stats(ii_stats-3).data,'ok')
            
            
        end
        
        
        
    end
    
    
end

for PACii=1:length(handles_out.drgbchoices.PACburstLowF)
    for percent_correct_ii=1:2
        figure((PACii-1)*2+percent_correct_ii)
        plot([-0.1 1],[-0.1 1],'-k')
        xlim([-0.1 1])
        ylim([-0.1 1])
        ylabel('AUC peak')
        xlabel('AUC trough')
        legend(handlesdrgb.drgbchoices.odorpair)
        title(['Area under the curve for Theta/' handles_out.drgbchoices.PACnames{PACii} ' ' handles_out.drgbchoices.group_no_names{groupNo}  ' ' handles_out.drgbchoices.per_lab{percent_correct_ii}])
    end
end

%Now do dimensonality

for figNo=7:12
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    set(hFig, 'units','normalized','position',[.2 .3 .6 .3])
    hold on
end



for fileNo=1:handlesdrgb.drgbchoices.no_files
    
    pname=handlesdrgb.drgbchoices.PathName{fileNo};
    fname=handlesdrgb.drgbchoices.FileName{fileNo};
    
    discriminant_name=[pname fname];
    load(discriminant_name)
    
    %Plot average percent correct for the LDA for peak and trough for
    %wavelet power referenced to PAC phase
    t=handles_out.t_power;
    
    
    
    for PACii=1:length(handles_out.drgbchoices.PACburstLowF)
        p_correct_stats=[];
        ii_stats=0;
        p_dim_stats=[];
        ii_dim_stats=0;
        glm_ii=0;
        glm_correct=[];
        glm_dim_ii=0;
        glm_dim=[];
        for percent_correct_ii=1:2
            
            bar_no=1+8*(fileNo-1);
            
            %Gather all the data
            no_mice=0;
            no_mice_included=0;

            all_dimensionality_peak=zeros(max(handles_out.drgbchoices.mouse_no),length(t));
            all_dimensionality_trough=zeros(max(handles_out.drgbchoices.mouse_no),length(t));
            
            for mouseNo=1:length(handles_out.discriminant_PACwavepower)
                try
                    if handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated==1
                        per_ii=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).no_trials;
                        no_mice=no_mice+1;
                        if (per_ii>=20)&(sum(no_mice==mice_excluded)==0)
                            no_mice_included=no_mice_included+1;
                            all_dimensionality_peak(no_mice_included,:)=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).dimensionality_peak;
                            all_dimensionality_trough(no_mice_included,:)=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).dimensionality_trough;
                        end
                    end
                catch
                end
            end
            
            
            no_mice_per(percent_correct_ii)=no_mice_included;
            %Plot dimensionality and save the data for
            %the ranksum
            
            
            %Now plot the pre-odor dimensionality and odor-induced change in dimensionality for the trough
            all_dimensionality_trough=all_dimensionality_trough(1:no_mice_included,:);
            
            data_pre=[];
            data_odor=[];
            for mouseNo=1:no_mice_included
                data_pre=[data_pre mean(all_dimensionality_trough(mouseNo,(t>=window_start(1))&(t<=window_end(1))),2)];
                data_odor=[data_odor mean(all_dimensionality_trough(mouseNo,(t>=window_start(2))&(t<=window_end(2))),2)];
            end
            
            %Plot pre
            figure(6+(PACii-1)*2+(percent_correct_ii-1)+1)
            hold on
            
            mean_dim_pre_trough=mean(data_pre',1)';
            CIpre_trough = bootci(1000, {@mean, data_pre})';
            h_thpre=bar(bar_no,mean_dim_pre_trough,'EdgeColor',[0.7 0.7 1],'FaceColor',[0.7 0.7 1.0]);
            plot([bar_no bar_no],CIpre_trough,'-k','LineWidth',3)
            
            bar_no=bar_no+1;
            
            %Plot odor
            mean_dim_odor_trough=mean(data_odor',1)';
            CIodor_trough = bootci(1000, {@mean, data_odor})';
            h_throd=bar(bar_no,mean_dim_odor_trough,'b');
            plot([bar_no bar_no],CIodor_trough,'-k','LineWidth',3)
            
            
            %                 ii_dim_stats=ii_dim_stats+1;
            %                 dim_stats(ii_dim_stats).data=data;
            %                 dim_stats(ii_dim_stats).description=[handles_out.drgbchoices.group_no_names{groupNo} ' ' ...
            %                     handles_out.drgbchoices.per_lab{percent_correct_ii}...
            %                     ' trough'];
            %                 dim_stats(ii_dim_stats).winNo=winNo;
            %                 dim_stats(ii_dim_stats).PACii=PACii;
            %                 dim_stats(ii_dim_stats).per_corr_ii=percent_correct_ii;
            %                 dim_stats(ii_dim_stats).groupNo=groupNo;
            %
            %                 glm_dim.data(glm_dim_ii+1:glm_dim_ii+length(data))=data;
            %                 glm_dim.PACii(glm_dim_ii+1:glm_dim_ii+length(data))=PACii;
            %                 glm_dim.group(glm_dim_ii+1:glm_dim_ii+length(data))=groupNo;
            %                 glm_dim.perCorr(glm_dim_ii+1:glm_dim_ii+length(data))=percent_correct_ii;
            %                 glm_dim.peak(glm_dim_ii+1:glm_dim_ii+length(data))=0;
            %                 glm_dim_ii=glm_dim_ii+length(data);
            
            
            %Now plot the percent correct for the peak
            all_dimensionality_peak=all_dimensionality_peak(1:no_mice_included,:);
            
            data_pre=[];
            data_odor=[];
            for mouseNo=1:no_mice_included
                data_pre=[data_pre mean(all_dimensionality_peak(mouseNo,(t>=window_start(1))&(t<=window_end(1))),2)];
                data_odor=[data_odor mean(all_dimensionality_peak(mouseNo,(t>=window_start(2))&(t<=window_end(2))),2)];
            end
            
            
            
            %Plot pre
            bar_no=bar_no+2;
            mean_dim_pre_peak=mean(data_pre',1)';
            CIpre_peak = bootci(1000, {@mean, data_pre})';
            h_pkpre=bar(bar_no,mean_dim_pre_peak,'EdgeColor',[1 0.7 0.7],'FaceColor',[1 0.7 0.7]);
            plot([bar_no bar_no],CIpre_peak,'-k','LineWidth',3)
            
            %Plot odor
            bar_no=bar_no+1;
            mean_dim_odor_peak=mean(data_odor',1)';
            CIodor_peak = bootci(1000, {@mean, data_odor})';
            h_pkod=bar(bar_no,mean_dim_odor_peak,'r');
            plot([bar_no bar_no],CIodor_peak,'-k','LineWidth',3)
            
            
            
            %                 data=[];
            %                 for mouseNo=1:no_mice_included
            %                     data=[data mean(all_dimensionality_peak(mouseNo,(t>=window_start(winNo))&(t<=window_end(winNo))),2)];
            %                 end
            %                 ii_dim_stats=ii_dim_stats+1;
            %                 dim_stats(ii_dim_stats).data=data;
            %                 dim_stats(ii_dim_stats).description=[handles_out.drgbchoices.group_no_names{groupNo} ' ' ...
            %                     handles_out.drgbchoices.per_lab{percent_correct_ii}...
            %                     ' peak'];
            %                 dim_stats(ii_dim_stats).winNo=winNo;
            %                 dim_stats(ii_dim_stats).PACii=PACii;
            %                 dim_stats(ii_dim_stats).per_corr_ii=percent_correct_ii;
            %                 dim_stats(ii_dim_stats).groupNo=groupNo;
            %
            %                 glm_dim.data(glm_dim_ii+1:glm_dim_ii+length(data))=data;
            %                 glm_dim.PACii(glm_dim_ii+1:glm_dim_ii+length(data))=PACii;
            %                 glm_dim.group(glm_dim_ii+1:glm_dim_ii+length(data))=groupNo;
            %                 glm_dim.perCorr(glm_dim_ii+1:glm_dim_ii+length(data))=percent_correct_ii;
            %                 glm_dim.peak(glm_dim_ii+1:glm_dim_ii+length(data))=1;
            %                 glm_dim_ii=glm_dim_ii+length(data);
            
            
            
            
        end
        
        
        
    end
    
   pfft=1;
   
end
 
for PACii=1:length(handles_out.drgbchoices.PACburstLowF)
    for percent_correct_ii=1:2
        %Odor on markers
        figure(6+(PACii-1)*2+(percent_correct_ii-1)+1)
        legend([h_thpre h_throd h_pkpre h_pkod],{'Trough pre-odor','Trough odor','Peak pre-odor','Peak odor'})
        xticks([3:8:43]) 
        xticklabels(handlesdrgb.drgbchoices.odorpair)
        title(['Dimensionality for Theta/' handles_out.drgbchoices.PACnames{PACii} ' ' handles_out.drgbchoices.group_no_names{groupNo}  ' ' handles_out.drgbchoices.per_lab{percent_correct_ii}])
        ylabel('Dimensionality')
    end
end
pffft=1;

function drgAnalyzeLFPDiscriminantBatchCaMKIIinf
%Analyzes the linear discriminant analysis performed by drgLFPDiscriminantBatch
%Takes as a in input the 'Discriminant_*.mat' output file from drgLFPDiscriminantBatch
%Performs an analysis of the timecourse for percent correct for LDA and for
%the PCA


warning('off')
close all
clear all


 
t_odor_arrival=0.1;

which_display=3;
mice_excluded=[];

[fname,pname,nCancel] = uigetfile({'Discriminant_*.mat'},'Select the perceptron LFP batch output file ...');
if nCancel
    inputPath = [pname,fname];
    pnameStart = pname;
    %save('timeCorr_cfg.mat','pnameStart','-append');
else
    error('Cancelled')
end

discriminant_name=[pname fname];
load(discriminant_name)

t_pac=handles_out.t_power';
t_power=handles_out.t_power';
handles.drg=handles_out.drg;

figNo=0;



% %Define the windows for analysis
% window_start=[-1 1.5];
% window_end=[0 2.5];
% no_wins=2;

%This is the window for odor
%
odor_from=1.5;
odor_to=2.5;


%This is the window for reinforcement
reinf_from=4.4;
reinf_to=5;

%This is the window for the AUC
AUC_from=0.5;
AUC_to=2.5;


%Plot average percent correct for the LDA for peak and trough for
%wavelet power referenced to PAC phase

%Note that I am hard setting this one
dt_tPRP=0.03333;
t=[-2:dt_tPRP:5];

% t=handles_out.t_power;
pcorr_out=[];
lick_out=[];

if length(handles_out.drgbchoices.PACnames)==3
    these_PACii=[1 3];
else
    these_PACii=[1 2];
end

for PACii=these_PACii
    
    p_correct_stats=[];
    ii_stats=0;
    p_corr_dt_peak_stats=[];
    ii_p_stats=0;
    p_corr_dt_trough_stats=[];
    ii_t_stats=0;
    glm_ii_pcdt=0;
    glm_pcorr_dt=[];

    glm_pcorr_dt_peak=[];
    glm_ii_pcdt_p=0;
    glm_pcorr_dt_trough=[];
    glm_ii_pcdt_t=0;
    
     
    
    for percent_correct_ii=1:2
        
        for groupNo=1:max(handles_out.drgbchoices.group_no)
            
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
                            pcorr_out.PACii(PACii).pcorr(percent_correct_ii).group(groupNo).mouseNos(no_mice_included)=mouseNo;
                            all_discriminant_correct_peak(no_mice_included,:)=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_peak;
                            all_discriminant_correct_trough(no_mice_included,:)=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_trough;
                            all_discriminant_correct_shuffled_peak(no_mice_included,:)=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_shuffled_peak;
                            all_discriminant_correct_shuffled_trough(no_mice_included,:)=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_shuffled_trough;
                            all_dimensionality_peak(no_mice_included,:)=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).dimensionality_peak;
                            all_dimensionality_trough(no_mice_included,:)=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).dimensionality_trough;
                        end
                    end
                catch
                end
            end
            
            fprintf(1, ['The number of mice included in the LDA analysis for this odor pair is %d\n\n\n'], no_mice_included)
            
            no_mice_per(percent_correct_ii)=no_mice_included;
            
            %Plot percent correct for the LDA and save the data for
            %the ranksum
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            hFig=figure(figNo);
            
            set(hFig, 'units','normalized','position',[.61 .4 .38 .38])
            
            hold on
            
            ax=gca;ax.LineWidth=3;
            
            %Note that I merge the shuffled for both peak and
            %trough for the average plot
            all_discriminant_correct_shuffled=zeros(2*no_mice_included,length(t));
            all_discriminant_correct_shuffled(1:no_mice_included,:)=all_discriminant_correct_shuffled_peak(1:no_mice_included,:);
            all_discriminant_correct_shuffled(no_mice_included+1:end,:)=all_discriminant_correct_shuffled_trough(1:no_mice_included,:);
            mean_dcsh=mean(all_discriminant_correct_shuffled,1)';
            if size(all_discriminant_correct_shuffled,1)>2
                CIdcsh = bootci(1000, {@mean, all_discriminant_correct_shuffled})';
                CIdcsh(:,1)=mean_dcsh-CIdcsh(:,1);
                CIdcsh(:,2)=CIdcsh(:,2)-mean_dcsh;
                [hlCR, hpCR] = boundedline(t,mean_dcsh, CIdcsh, 'k');
            else
                plot(t,mean_dcsh,'-k')
            end
            
            data=[];
            for mouseNo=1:no_mice_included
                data=[data mean(all_discriminant_correct_shuffled_peak(mouseNo,(t>=odor_from)&(t<=odor_to)),2)];
            end
            
            dataAUC=[];
            for mouseNo=1:no_mice_included
                dataAUC=[dataAUC (mean(all_discriminant_correct_shuffled_peak(mouseNo,(t>=AUC_from)&(t<=AUC_to)),2)-50)/50];
            end
     
            pcorr_out.PACii(PACii).shuffled_peak.pcorr(percent_correct_ii).group(groupNo).odor_data=data;
            
            glm_pcorr_dt.data(glm_ii_pcdt+1:glm_ii_pcdt+length(data))=data;
            glm_pcorr_dt.dataAUC(glm_ii_pcdt+1:glm_ii_pcdt+length(data))=dataAUC;
            glm_pcorr_dt.group(glm_ii_pcdt+1:glm_ii_pcdt+length(data))=groupNo;
            glm_pcorr_dt.perCorr(glm_ii_pcdt+1:glm_ii_pcdt+length(data))=percent_correct_ii;
            %             glm_pcorr_dt.shuffled(glm_ii_pcdt+1:glm_ii_pcdt+length(data))=1;
            glm_pcorr_dt.pts(glm_ii_pcdt+1:glm_ii_pcdt+length(data))=3;
            glm_ii_pcdt=glm_ii_pcdt+length(data);
            
      
            
            datat=[];
            for mouseNo=1:no_mice_included
                data=[data mean(all_discriminant_correct_shuffled_trough(mouseNo,(t>=odor_from)&(t<=odor_to)),2)];
                datat=[datat mean(all_discriminant_correct_shuffled_trough(mouseNo,(t>=odor_from)&(t<=odor_to)),2)];
            end
            
            for mouseNo=1:no_mice_included
                dataAUC=[dataAUC (mean(all_discriminant_correct_shuffled_trough(mouseNo,(t>=AUC_from)&(t<=AUC_to)),2)-50)/50];
            end
            
            pcorr_out.PACii(PACii).shuffled_trough.pcorr(percent_correct_ii).group(groupNo).odor_data=datat;
            pcorr_out.PACii(PACii).shuffled.pcorr(percent_correct_ii).group(groupNo).odor_data=data;
            
            ii_stats=ii_stats+1;
            p_correct_stats(ii_stats).data=data;
            p_correct_stats(ii_stats).description=[handles_out.drgbchoices.group_no_names{groupNo}  ' ' handles_out.drgbchoices.per_lab{percent_correct_ii} ' shuffled'];
            
            
            glm_pcorr_dt.data(glm_ii_pcdt+1:glm_ii_pcdt+length(data))=data;
            glm_pcorr_dt.dataAUC(glm_ii_pcdt+1:glm_ii_pcdt+length(data))=dataAUC;
            glm_pcorr_dt.group(glm_ii_pcdt+1:glm_ii_pcdt+length(data))=groupNo;
            glm_pcorr_dt.perCorr(glm_ii_pcdt+1:glm_ii_pcdt+length(data))=percent_correct_ii;
            %             glm_pcorr_dt.shuffled(glm_ii_pcdt+1:glm_ii_pcdt+length(data))=1;
            glm_pcorr_dt.pts(glm_ii_pcdt+1:glm_ii_pcdt+length(data))=3;
            glm_ii_pcdt=glm_ii_pcdt+length(data);
            
  
            
            
            %Now plot the percent correct for the trough
            all_discriminant_correct_trough=all_discriminant_correct_trough(1:no_mice_included,:);
            mean_dc_trough=mean(all_discriminant_correct_trough,1)';
            if size(all_discriminant_correct_trough,1)>2
                CIdc_trough = bootci(1000, {@mean, all_discriminant_correct_trough})';
                CIdc_trough(:,1)=mean_dc_trough-CIdc_trough(:,1);
                CIdc_trough(:,2)=CIdc_trough(:,2)-mean_dc_trough;
                [hlCR, hpCR] = boundedline(t,mean_dc_trough, CIdc_trough, 'b');
            else
                plot(t,mean_dc_trough,'b')
            end
            
            
            data=[];
            for mouseNo=1:no_mice_included
                data=[data mean(all_discriminant_correct_trough(mouseNo,(t>=odor_from)&(t<=odor_to)),2)];
            end
            
             dataAUC=[];
            for mouseNo=1:no_mice_included
                dataAUC=[dataAUC (mean(all_discriminant_correct_trough(mouseNo,(t>=AUC_from)&(t<=AUC_to)),2)-50)/50];
            end
            
            pcorr_out.PACii(PACii).trough.pcorr(percent_correct_ii).group(groupNo).odor_data=data;
            pcorr_out.PACii(PACii).trough.pcorr(percent_correct_ii).group(groupNo).odor_dataAUC=dataAUC;
            
            ii_stats=ii_stats+1;
            p_correct_stats(ii_stats).data=data;
            p_correct_stats(ii_stats).description=[handles_out.drgbchoices.group_no_names{groupNo} ' ' handles_out.drgbchoices.per_lab{percent_correct_ii} ' trough'];
            
            
            glm_pcorr_dt.data(glm_ii_pcdt+1:glm_ii_pcdt+length(data))=data;
            glm_pcorr_dt.dataAUC(glm_ii_pcdt+1:glm_ii_pcdt+length(data))=dataAUC;
            glm_pcorr_dt.group(glm_ii_pcdt+1:glm_ii_pcdt+length(data))=groupNo;
            glm_pcorr_dt.perCorr(glm_ii_pcdt+1:glm_ii_pcdt+length(data))=percent_correct_ii;
            %             glm_pcorr_dt.shuffled(glm_ii_pcdt+1:glm_ii_pcdt+length(data))=0;
            glm_pcorr_dt.pts(glm_ii_pcdt+1:glm_ii_pcdt+length(data))=2;
            glm_ii_pcdt=glm_ii_pcdt+length(data);
            

            
            glm_pcorr_dt_trough.data(glm_ii_pcdt_t+1:glm_ii_pcdt_t+length(data))=data;
            glm_pcorr_dt_trough.dataAUC(glm_ii_pcdt_t+1:glm_ii_pcdt_t+length(data))=dataAUC;
            glm_pcorr_dt_trough.group(glm_ii_pcdt_t+1:glm_ii_pcdt_t+length(data))=groupNo;
            glm_pcorr_dt_trough.perCorr(glm_ii_pcdt_t+1:glm_ii_pcdt_t+length(data))=percent_correct_ii;
            glm_ii_pcdt_t=glm_ii_pcdt_t+length(data);
            
            ii_t_stats=ii_t_stats+1;
            p_corr_dt_trough_stats(ii_t_stats).data=data;
            p_corr_dt_trough_stats(ii_t_stats).description=[handles_out.drgbchoices.group_no_names{groupNo} ' ' ...
                handles_out.drgbchoices.per_lab{percent_correct_ii}];
            
            
            data=[];
            for mouseNo=1:no_mice_included
                data=[data mean(all_discriminant_correct_trough(mouseNo,(t>=reinf_from)&(t<=reinf_to)),2)];
            end
            
            pcorr_out.PACii(PACii).trough.pcorr(percent_correct_ii).group(groupNo).reinf_data=data;

            
            
            %Now plot the percent correct for the peak
            all_discriminant_correct_peak=all_discriminant_correct_peak(1:no_mice_included,:);
            mean_dc_peak=mean(all_discriminant_correct_peak,1)';
            if size(all_discriminant_correct_peak,1)>2
                CIdc_peak = bootci(1000, {@mean, all_discriminant_correct_peak})';
                CIdc_peak(:,1)=mean_dc_peak-CIdc_peak(:,1);
                CIdc_peak(:,2)=CIdc_peak(:,2)-mean_dc_peak;
                [hlCR, hpCR] = boundedline(t,mean_dc_peak, CIdc_peak, 'r');
            else
                plot(t,mean_dc_peak,'r')
            end
            
            
            data=[];
            for mouseNo=1:no_mice_included
                data=[data mean(all_discriminant_correct_peak(mouseNo,(t>=odor_from)&(t<=odor_to)),2)];
            end
            
            dataAUC=[];
            for mouseNo=1:no_mice_included
                dataAUC=[dataAUC (mean(all_discriminant_correct_peak(mouseNo,(t>=AUC_from)&(t<=AUC_to)),2)-50)/50];
            end
            
            pcorr_out.PACii(PACii).peaks.pcorr(percent_correct_ii).group(groupNo).odor_data=data;
            pcorr_out.PACii(PACii).peaks.pcorr(percent_correct_ii).group(groupNo).odor_dataAUC=dataAUC;
            
            ii_stats=ii_stats+1;
            p_correct_stats(ii_stats).data=data;
            p_correct_stats(ii_stats).description=[handles_out.drgbchoices.group_no_names{groupNo} ' ' ...
                handles_out.drgbchoices.per_lab{percent_correct_ii}...
                ' peak'];
            
            
            glm_pcorr_dt.data(glm_ii_pcdt+1:glm_ii_pcdt+length(data))=data;
            glm_pcorr_dt.dataAUC(glm_ii_pcdt+1:glm_ii_pcdt+length(data))=dataAUC;
            glm_pcorr_dt.PACii(glm_ii_pcdt+1:glm_ii_pcdt+length(data))=PACii;
            glm_pcorr_dt.group(glm_ii_pcdt+1:glm_ii_pcdt+length(data))=groupNo;
            glm_pcorr_dt.perCorr(glm_ii_pcdt+1:glm_ii_pcdt+length(data))=percent_correct_ii;
            %             glm_pcorr_dt.shuffled(glm_ii_pcdt+1:glm_ii_pcdt+length(data))=0;
            glm_pcorr_dt.pts(glm_ii_pcdt+1:glm_ii_pcdt+length(data))=1;
            glm_ii_pcdt=glm_ii_pcdt+length(data);
            

            
            glm_pcorr_dt_peak.data(glm_ii_pcdt_p+1:glm_ii_pcdt_p+length(data))=data;
            glm_pcorr_dt_peak.dataAUC(glm_ii_pcdt_p+1:glm_ii_pcdt_p+length(data))=dataAUC;
            glm_pcorr_dt_peak.group(glm_ii_pcdt_p+1:glm_ii_pcdt_p+length(data))=groupNo;
            glm_pcorr_dt_peak.perCorr(glm_ii_pcdt_p+1:glm_ii_pcdt_p+length(data))=percent_correct_ii;
            glm_ii_pcdt_p=glm_ii_pcdt_p+length(data);
            
            ii_p_stats=ii_p_stats+1;
            p_corr_dt_peak_stats(ii_p_stats).data=data;
            p_corr_dt_peak_stats(ii_p_stats).description=[handles_out.drgbchoices.group_no_names{groupNo} ' ' ...
                handles_out.drgbchoices.per_lab{percent_correct_ii}];
            
            data=[];
            for mouseNo=1:no_mice_included
                data=[data mean(all_discriminant_correct_peak(mouseNo,(t>=reinf_from)&(t<=reinf_to)),2)];
            end
            
            pcorr_out.PACii(PACii).peaks.pcorr(percent_correct_ii).group(groupNo).reinf_data=data;
            
            %Odor on markers
            plot([0 0],[0 100],'-k')
            odorhl=plot([0 2.5],[10 10],'-k','LineWidth',5);
            plot([2.5 2.5],[0 100],'-k')
            
            title(['LDA for theta/' handles_out.drgbchoices.PACnames{PACii} ' ' handles_out.drgbchoices.group_no_names{groupNo}  ' ' handles_out.drgbchoices.per_lab{percent_correct_ii}])
            
            xlabel('Time (sec)')
            ylabel(['% correct '  handles_out.drgbchoices.per_lab{percent_correct_ii}])
            legend('Shuffled','Trough','Peak')
           
            
%             if figNo==13
%                pfft=0; 
%             end
            
        end
        
        
    end
    
    %Bar graph plot for odor for peak percent correct
    %Plot the average
    figNo = figNo +1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    
    set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
    hold on
    
    ax=gca;ax.LineWidth=3;
    
    bar_offset = 0;
    
    for per_ii=2:-1:1
        
        for grNo=1:max(handles_out.drgbchoices.group_no)
            bar_offset = bar_offset +1;
            
            
            switch grNo
                case 1
                    bar(bar_offset,mean(pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data),'g','LineWidth', 3,'EdgeColor','none')
                case 2
                    bar(bar_offset,mean(pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data),'b','LineWidth', 3,'EdgeColor','none')
                case 3
                    bar(bar_offset,mean(pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data),'y','LineWidth', 3,'EdgeColor','none')
            end
            
            try
                CI = bootci(1000, {@mean, pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data},'type','cper');
                plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                plot(bar_offset*ones(1,length(pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data)),pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
            catch
            end
            
        end
        bar_offset = bar_offset + 2;
        
    end
    
    ylim([45 100])
    
    title(['LDA percent correct for odor for peak theta/' handles_out.drgbchoices.PACnames{PACii} ])
    
    
    xticks([1 2 3 5 6 7])
    xticklabels({'nwt', 'nH', 'nKO', 'pwt', 'pH', 'pKO'})
    
    
    ylabel('Percent correct')
    
    %Bar graph plot for odor for peak percent correct AUC
    %Plot the average
    figNo = figNo +1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    
    set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
    hold on
    
    ax=gca;ax.LineWidth=3;
    
    bar_offset = 0;
    
    for per_ii=2:-1:1
        
        for grNo=1:max(handles_out.drgbchoices.group_no)
            bar_offset = bar_offset +1;
            
            
            switch grNo
                case 1
                    bar(bar_offset,mean(pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_dataAUC),'g','LineWidth', 3,'EdgeColor','none')
                case 2
                    bar(bar_offset,mean(pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_dataAUC),'b','LineWidth', 3,'EdgeColor','none')
                case 3
                    bar(bar_offset,mean(pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_dataAUC),'y','LineWidth', 3,'EdgeColor','none')
            end
            
            try
                CI = bootci(1000, {@mean, pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_dataAUC},'type','cper');
                plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                plot(bar_offset*ones(1,length(pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_dataAUC)),pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_dataAUC,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
            catch
            end
            
        end
        bar_offset = bar_offset + 2;
        
    end
    
    ylim([0 0.9])
    
    title(['LDA AUC percent correct for odor for peak theta/' handles_out.drgbchoices.PACnames{PACii} ])
    
    
    xticks([1 2 3 5 6 7])
    xticklabels({'nwt', 'nH', 'nKO', 'pwt', 'pH', 'pKO'})
    
    
    ylabel('Percent correct')
    
    %Bar graph plot for reinforcement for peak percent correct
    %Plot the average
    figNo = figNo +1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    
    set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
    hold on
    
    ax=gca;ax.LineWidth=3;
    
    bar_offset = 0;
    
    for per_ii=2:-1:1
        
        for grNo=1:max(handles_out.drgbchoices.group_no)
            bar_offset = bar_offset +1;
            
            
            switch grNo
                case 1
                    bar(bar_offset,mean(pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).reinf_data),'g','LineWidth', 3,'EdgeColor','none')
                case 2
                    bar(bar_offset,mean(pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).reinf_data),'b','LineWidth', 3,'EdgeColor','none')
                case 3
                    bar(bar_offset,mean(pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).reinf_data),'y','LineWidth', 3,'EdgeColor','none')
            end
            
            try
                CI = bootci(1000, {@mean, pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).reinf_data},'type','cper');
                plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                plot(bar_offset*ones(1,length(pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).reinf_data)),pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).reinf_data,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
            catch
            end
            
        end
        bar_offset = bar_offset + 2;
        
    end
    
    ylim([45 100])
    
    title(['LDA percent correct for reinforcement for peak theta/' handles_out.drgbchoices.PACnames{PACii} ])
    
    
    xticks([1 2 3 5 6 7])
    xticklabels({'nwt', 'nH', 'nKO', 'pwt', 'pH', 'pKO'})
    
    
    ylabel('Percent correct')
    
    %Bar graph plot for odor for trough percent correct
    %Plot the average
    figNo = figNo +1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    
    set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
    hold on
    
    ax=gca;ax.LineWidth=3;
    
    bar_offset = 0;
    
    for per_ii=2:-1:1
        
        for grNo=1:max(handles_out.drgbchoices.group_no)
            bar_offset = bar_offset +1;
            
            
            switch grNo
                case 1
                    bar(bar_offset,mean(pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data),'g','LineWidth', 3,'EdgeColor','none')
                case 2
                    bar(bar_offset,mean(pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data),'b','LineWidth', 3,'EdgeColor','none')
                case 3
                    bar(bar_offset,mean(pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data),'y','LineWidth', 3,'EdgeColor','none')
            end
            
            try
                CI = bootci(1000, {@mean, pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data},'type','cper');
                plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                plot(bar_offset*ones(1,length(pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data)),pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
            catch
            end
            
        end
        bar_offset = bar_offset + 2;
        
    end
    
    ylim([45 100])
    
    title(['LDA percent correct for odor for trough theta/' handles_out.drgbchoices.PACnames{PACii} ])
    
    
    xticks([1 2 3 5 6 7])
    xticklabels({'nwt', 'nH', 'nKO', 'pwt', 'pH', 'pKO'})
    
    
    ylabel('Percent correct')
    
    %Bar graph plot for odor for trough percent correct AUC
    %Plot the average
    figNo = figNo +1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    
    set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
    hold on
    
    ax=gca;ax.LineWidth=3;
    
    bar_offset = 0;
    
    for per_ii=2:-1:1
        
        for grNo=1:max(handles_out.drgbchoices.group_no)
            bar_offset = bar_offset +1;
            
            
            switch grNo
                case 1
                    bar(bar_offset,mean(pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_dataAUC),'g','LineWidth', 3,'EdgeColor','none')
                case 2
                    bar(bar_offset,mean(pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_dataAUC),'b','LineWidth', 3,'EdgeColor','none')
                case 3
                    bar(bar_offset,mean(pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_dataAUC),'y','LineWidth', 3,'EdgeColor','none')
            end
            
            try
                CI = bootci(1000, {@mean, pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_dataAUC},'type','cper');
                plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                plot(bar_offset*ones(1,length(pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_dataAUC)),pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_dataAUC,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
            catch
            end
            
        end
        bar_offset = bar_offset + 2;
        
    end
    
    ylim([45 100])
    
    title(['LDA AUC percent correct for odor for trough theta/' handles_out.drgbchoices.PACnames{PACii} ])
    
    
    xticks([1 2 3 5 6 7])
    xticklabels({'nwt', 'nH', 'nKO', 'pwt', 'pH', 'pKO'})
    
    
    ylabel('Percent correct')
    
     %Bar graph plot for reinforcment for trough percent correct
    %Plot the average
    figNo = figNo +1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    
    set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
    hold on
    
    ax=gca;ax.LineWidth=3;
    
    bar_offset = 0;
    
    for per_ii=2:-1:1
        
        for grNo=1:max(handles_out.drgbchoices.group_no)
            bar_offset = bar_offset +1;
            
            
            switch grNo
                case 1
                    bar(bar_offset,mean(pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).reinf_data),'g','LineWidth', 3,'EdgeColor','none')
                case 2
                    bar(bar_offset,mean(pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).reinf_data),'b','LineWidth', 3,'EdgeColor','none')
                case 3
                    bar(bar_offset,mean(pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).reinf_data),'y','LineWidth', 3,'EdgeColor','none')
            end
            
            try
                CI = bootci(1000, {@mean, pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).reinf_data},'type','cper');
                plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                plot(bar_offset*ones(1,length(pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).reinf_data)),pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).reinf_data,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
            catch
            end
            
        end
        bar_offset = bar_offset + 2;
        
    end
    
    ylim([45 100])
    
    title(['LDA percent correct for reinforcement for trough theta/' handles_out.drgbchoices.PACnames{PACii} ])
    
    
    xticks([1 2 3 5 6 7])
    xticklabels({'nwt', 'nH', 'nKO', 'pwt', 'pH', 'pKO'})
    
    
    ylabel('Percent correct')
    
    
    %Perform the glm
    fprintf(1, ['\n\nglm for mean percent correct peak Theta/' handles_out.drgbchoices.PACnames{PACii} '\n'])
    tbl = table(glm_pcorr_dt_peak.data',glm_pcorr_dt_peak.group',glm_pcorr_dt_peak.perCorr',...
        'VariableNames',{'mean_pc','group','naive_proficient'});
    mdl = fitglm(tbl,'mean_pc~group+naive_proficient+group*naive_proficient'...
        ,'CategoricalVars',[2,3])
    
    %Do ranksum/t test
    fprintf(1, ['\n\nRanksum or t-test p values for mean percent correct peak Theta/' handles_out.drgbchoices.PACnames{PACii} '\n'])
    try
        [output_data] = drgMutiRanksumorTtest(p_corr_dt_peak_stats);
        fprintf(1, '\n\n')
    catch
    end
    
    %Perform the glm
    fprintf(1, ['\n\nglm for mean percent correct trough Theta/' handles_out.drgbchoices.PACnames{PACii} '\n'])
    tbl = table(glm_pcorr_dt_trough.data',glm_pcorr_dt_trough.group',glm_pcorr_dt_trough.perCorr',...
        'VariableNames',{'mean_pc','group','naive_proficient'});
    mdl = fitglm(tbl,'mean_pc~group+naive_proficient+group*naive_proficient'...
        ,'CategoricalVars',[2,3])
    
    %Do ranksum/t test
    fprintf(1, ['\n\nRanksum or t-test p values for mean percent correct trough Theta/' handles_out.drgbchoices.PACnames{PACii} '\n'])
    try
        [output_data] = drgMutiRanksumorTtest(p_corr_dt_peak_stats);
        fprintf(1, '\n\n')
    catch
    end
    
    %Perform the glm
    fprintf(1, ['\n\nglm for mean percent correct for Theta/' handles_out.drgbchoices.PACnames{PACii} '\n'])
    tbl = table(glm_pcorr_dt.data',glm_pcorr_dt.group',glm_pcorr_dt.perCorr',glm_pcorr_dt.pts',...
        'VariableNames',{'mean_pc','group','perCorr','peak_trough_shuffled'});
    mdl = fitglm(tbl,'mean_pc~group+perCorr+peak_trough_shuffled+group*peak_trough_shuffled*perCorr'...
        ,'CategoricalVars',[2,3,4])
    
    %Do ranksum/t test
    fprintf(1, ['\n\nRanksum or t-test p values for mean percent correct  for Theta/' handles_out.drgbchoices.PACnames{PACii} '\n'])
    try
        [output_data] = drgMutiRanksumorTtest(p_correct_stats);
        fprintf(1, '\n\n')
    catch
    end
    
     %Perform the glm
    fprintf(1, ['\n\nglm for mean AUC percent correct for Theta/' handles_out.drgbchoices.PACnames{PACii} '\n'])
    tbl = table(glm_pcorr_dt.dataAUC',glm_pcorr_dt.group',glm_pcorr_dt.perCorr',glm_pcorr_dt.pts',...
        'VariableNames',{'mean_pc','group','perCorr','peak_trough_shuffled'});
    mdl = fitglm(tbl,'mean_pc~group+perCorr+peak_trough_shuffled+group*peak_trough_shuffled*perCorr'...
        ,'CategoricalVars',[2,3,4])
    
    
    
end

 
%Now plot log(p) and find decision times
disc_time=[];
for PACii=these_PACii
    
    
    glm_d_time=[];
    glm_d_t_ii=0;
    t_det_stats=[];
    ii_stats=0;
    
    
    
    for percent_correct_ii=1:2
        
        for groupNo=1:max(handles_out.drgbchoices.group_no)
            
            %Gather all the data
            no_mice=0;
            no_mice_included=0;
            all_discriminant_p_val_lick=zeros(max(handles_out.drgbchoices.mouse_no),length(t));
            all_discriminant_p_val_peak=zeros(max(handles_out.drgbchoices.mouse_no),length(t));
            all_discriminant_p_val_trough=zeros(max(handles_out.drgbchoices.mouse_no),length(t));
            
            for mouseNo=1:length(handles_out.discriminant_PACwavepower)
                try
                    if handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated==1
                        per_ii=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).no_trials;
                        no_mice=no_mice+1;
                        if (per_ii>=20)&(sum(no_mice==mice_excluded)==0)
                            no_mice_included=no_mice_included+1;
                            pval_out.PACii(PACii).pcorr(percent_correct_ii).group(groupNo).mouseNos(no_mice_included)=mouseNo;
                            disc_time.PACii(PACii).pcorr(percent_correct_ii).group(groupNo).mouseNos(no_mice_included)=mouseNo;
                            all_discriminant_p_val_lick(no_mice_included,:)=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).p_val_lick;
                            all_discriminant_p_val_peak(no_mice_included,:)=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).p_val_peak;
                            all_discriminant_p_val_trough(no_mice_included,:)=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).p_val_trough;
                        end
                    end
                catch
                end
            end
            
            
            
            no_mice_per(percent_correct_ii)=no_mice_included;
            
            %Plot the p value timecourse and compute the decision time
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            hFig=figure(figNo);
            
            set(hFig, 'units','normalized','position',[.61 .4 .38 .38])
            
            hold on
            
            ax=gca;ax.LineWidth=3;
            
            all_discriminant_p_val_lick=log10(all_discriminant_p_val_lick(1:no_mice_included,:));
            all_discriminant_p_val_peak=log10(all_discriminant_p_val_peak(1:no_mice_included,:));
            all_discriminant_p_val_trough=log10(all_discriminant_p_val_trough(1:no_mice_included,:));
            
            mean_p_val_licks=nanmean(all_discriminant_p_val_lick,1)';
            if size(all_discriminant_p_val_lick,1)>2
                CIlicks = bootci(1000, {@nanmean, all_discriminant_p_val_lick})';
                CIlicks(:,1)=mean_p_val_licks-CIlicks(:,1);
                CIlicks(:,2)=CIlicks(:,2)-mean_p_val_licks;
                [hlick, hpCR] = boundedline(t,mean_p_val_licks, CIlicks, 'k');
            else
                plot(t,mean_p_val_licks,'-k')
            end
            
            data=[];
            for mouseNo=1:no_mice_included
                data=[data mean(all_discriminant_p_val_lick(mouseNo,(t>=odor_from)&(t<=odor_to)),2)];
            end
            
            pval_out.PACii(PACii).lick.pcorr(percent_correct_ii).group(groupNo).odor_data=data;
            
            wing=1;
            t_detect=zeros(1,no_mice_included);
            jj_start=find(t>=t_odor_arrival,1,'first');
            %Find the discrimination time
            for ii=1:no_mice_included
                this_p_val_licks=[];
                this_p_val_licks=all_discriminant_p_val_lick(ii,:);
                found_disc_t=0;
                while (found_disc_t==0)&(jj_start<length(t))
                    ii_next=find(this_p_val_licks(1,jj_start:end)<=log10(0.05),1,'first');
                    if isempty(ii_next)
                        jj_start=length(t);
                        ii_next=1;
                    else
                        if jj_start+ii_next+wing>length(t)
                            found_disc_t=1;
                        else
                            if sum(this_p_val_licks(1,jj_start+ii_next:jj_start+ii_next+wing)>log10(0.05))>=1
                                jj_start=jj_start+ii_next;
                            else
                                found_disc_t=1;
                            end
                        end
                    end
                end
                t_detect(ii)=t(jj_start+ii_next-1)-t_odor_arrival;
            end
            
            fprintf(1, ['Lick discrimination time (sec) %d\n'],mean(t_detect))
            
            lick_t_detect=t_detect;
            
            ii_stats=ii_stats+1;
            t_det_stats(ii_stats).data=t_detect;
            t_det_stats(ii_stats).description=[handles_out.drgbchoices.group_no_names{groupNo} ' ' ...
                handles_out.drgbchoices.per_lab{percent_correct_ii}...
                ' licks'];
            
            glm_d_time.data(glm_d_t_ii+1:glm_d_t_ii+length(t_detect))=t_detect;
            glm_d_time.group(glm_d_t_ii+1:glm_d_t_ii+length(t_detect))=groupNo;
            glm_d_time.perCorr(glm_d_t_ii+1:glm_d_t_ii+length(t_detect))=percent_correct_ii;
            %             glm_pcorr_dt.shuffled(glm_ii_pcdt+1:glm_ii_pcdt+length(data))=0;
            glm_d_time.ptl(glm_d_t_ii+1:glm_d_t_ii+length(t_detect))=1;
            glm_d_t_ii=glm_d_t_ii+length(t_detect);
            
            disc_time.PACii(PACii).licks.pcorr(percent_correct_ii).group(groupNo).data=t_detect;
            
            mean_p_val_troughs=nanmean(all_discriminant_p_val_trough,1)';
            if size(all_discriminant_p_val_trough,1)>2
                CItroughs = bootci(1000, {@nanmean, all_discriminant_p_val_trough})';
                CItroughs(:,1)=mean_p_val_troughs-CItroughs(:,1);
                CItroughs(:,2)=CItroughs(:,2)-mean_p_val_troughs;
                [hltrough, hpCR] = boundedline(t,mean_p_val_troughs, CItroughs, 'b');
            else
                plot(t,mean_p_val_troughs,'-b')
            end
            
            data=[];
            for mouseNo=1:no_mice_included
                data=[data mean(all_discriminant_p_val_trough(mouseNo,(t>=odor_from)&(t<=odor_to)),2)];
            end
            
            pval_out.PACii(PACii).trough.pcorr(percent_correct_ii).group(groupNo).odor_data=data;
            
            t_detect=zeros(1,no_mice_included);
            jj_start=find(t>=t_odor_arrival,1,'first');
            %Find the discrimination time
            for ii=1:no_mice_included
                this_p_val_troughs=[];
                this_p_val_troughs=all_discriminant_p_val_trough(ii,:);
                found_disc_t=0;
                while (found_disc_t==0)&(jj_start<length(t))
                    ii_next=find(this_p_val_troughs(1,jj_start:end)<=log10(0.05),1,'first');
                    if isempty(ii_next)
                        jj_start=length(t);
                        ii_next=1;
                    else
                        if jj_start+ii_next+wing>length(t)
                            found_disc_t=1;
                        else
                            if sum(this_p_val_troughs(1,jj_start+ii_next:jj_start+ii_next+wing)>log10(0.05))>=1
                                jj_start=jj_start+ii_next;
                            else
                                found_disc_t=1;
                            end
                        end
                    end
                end
                t_detect(ii)=t(jj_start+ii_next-1)-t_odor_arrival;
            end
            
            trough_t_detect=t_detect;
            
            fprintf(1, ['Trough discrimination time (sec) %d\n'],mean(t_detect))
            
            ii_stats=ii_stats+1;
            t_det_stats(ii_stats).data=t_detect;
            t_det_stats(ii_stats).description=[handles_out.drgbchoices.group_no_names{groupNo} ' ' ...
                handles_out.drgbchoices.per_lab{percent_correct_ii}...
                ' troughs'];
            
            disc_time.PACii(PACii).trough.pcorr(percent_correct_ii).group(groupNo).data=t_detect;
            
            glm_d_time.data(glm_d_t_ii+1:glm_d_t_ii+length(t_detect))=t_detect;
            glm_d_time.group(glm_d_t_ii+1:glm_d_t_ii+length(t_detect))=groupNo;
            glm_d_time.perCorr(glm_d_t_ii+1:glm_d_t_ii+length(t_detect))=percent_correct_ii;
            %             glm_pcorr_dt.shuffled(glm_ii_pcdt+1:glm_ii_pcdt+length(data))=0;
            glm_d_time.ptl(glm_d_t_ii+1:glm_d_t_ii+length(t_detect))=2;
            glm_d_t_ii=glm_d_t_ii+length(t_detect);
            
            mean_p_val_peaks=nanmean(all_discriminant_p_val_peak,1)';
            if size(all_discriminant_p_val_peak,1)>2
                CIpeaks = bootci(1000, {@nanmean, all_discriminant_p_val_peak})';
                CIpeaks(:,1)=mean_p_val_peaks-CIpeaks(:,1);
                CIpeaks(:,2)=CIpeaks(:,2)-mean_p_val_peaks;
                [hpeak, hpCR] = boundedline(t,mean_p_val_peaks, CIpeaks, 'r');
            else
                plot(t,mean_p_val_peaks,'-r')
            end
            
            data=[];
            for mouseNo=1:no_mice_included
                data=[data mean(all_discriminant_p_val_peak(mouseNo,(t>=odor_from)&(t<=odor_to)),2)];
            end
            
            pval_out.PACii(PACii).peaks.pcorr(percent_correct_ii).group(groupNo).odor_data=data;
            
            t_detect=zeros(1,no_mice_included);
            jj_start=find(t>=t_odor_arrival,1,'first');
            %Find the discrimination time
            for ii=1:no_mice_included
                this_p_val_peaks=[];
                this_p_val_peaks=all_discriminant_p_val_peak(ii,:);
                found_disc_t=0;
                while (found_disc_t==0)&(jj_start<length(t))
                    ii_next=find(this_p_val_peaks(1,jj_start:end)<=log10(0.05),1,'first');
                    if isempty(ii_next)
                        jj_start=length(t);
                        ii_next=1;
                    else
                        if jj_start+ii_next+wing>length(t)
                            found_disc_t=1;
                        else
                            if sum(this_p_val_peaks(1,jj_start+ii_next:jj_start+ii_next+wing)>log10(0.05))>=1
                                jj_start=jj_start+ii_next;
                            else
                                found_disc_t=1;
                            end
                        end
                    end
                end
                t_detect(ii)=t(jj_start+ii_next-1)-t_odor_arrival;
            end
            
            fprintf(1, ['Peak discrimination time (sec) %d\n'],mean(t_detect))
            
            peak_t_detect=t_detect;
            
            ii_stats=ii_stats+1;
            t_det_stats(ii_stats).data=t_detect;
            t_det_stats(ii_stats).description=[handles_out.drgbchoices.group_no_names{groupNo} ' ' ...
                handles_out.drgbchoices.per_lab{percent_correct_ii}...
                ' peaks'];
            disc_time.PACii(PACii).peaks.pcorr(percent_correct_ii).group(groupNo).data=t_detect;
            
            glm_d_time.data(glm_d_t_ii+1:glm_d_t_ii+length(t_detect))=t_detect;
            glm_d_time.group(glm_d_t_ii+1:glm_d_t_ii+length(t_detect))=groupNo;
            glm_d_time.perCorr(glm_d_t_ii+1:glm_d_t_ii+length(t_detect))=percent_correct_ii;
            %             glm_pcorr_dt.shuffled(glm_ii_pcdt+1:glm_ii_pcdt+length(data))=0;
            glm_d_time.ptl(glm_d_t_ii+1:glm_d_t_ii+length(t_detect))=3;
            glm_d_t_ii=glm_d_t_ii+length(t_detect);
            
            plot(t,mean_p_val_licks,'-k')
            plot(t,mean_p_val_troughs,'-b')
            plot(t,mean_p_val_peaks,'-r')
            plot([t(1) t(end)],[log10(0.05) log10(0.05)],'-r')
            
            
            %Odor on markers
            
            odorhl=plot([0 2.5],[0 0],'-k','LineWidth',5);
            
            ylim([-80 0])
            title(['p value for Theta/' handles_out.drgbchoices.PACnames{PACii} ' ' handles_out.drgbchoices.group_no_names{groupNo}  ' ' handles_out.drgbchoices.per_lab{percent_correct_ii}])
            
            xlabel('Time (sec)')
            ylabel(['p value '  handles_out.drgbchoices.per_lab{percent_correct_ii}])
            
            try
                legend([hlick hltrough hpeak],{'Licks','Trough','Peak'})
            catch
            end
            
            
            
        end
    end
    
    if PACii==1
        %Bar graph plot for lick discrimination time
        %Plot the average
        figNo = figNo +1;
        
        try
            close(figNo)
        catch
        end
        hFig=figure(figNo);
        
        
        set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
        hold on
        
        ax=gca;ax.LineWidth=3;
        
        bar_offset = 0;
        
        for per_ii=2:-1:1
            
            for grNo=1:max(handles_out.drgbchoices.group_no)
                bar_offset = bar_offset +1;
                
                
                switch grNo
                    case 1
                        bar(bar_offset,mean(disc_time.PACii(PACii).licks.pcorr(per_ii).group(grNo).data),'g','LineWidth', 3,'EdgeColor','none')
                    case 2
                        bar(bar_offset,mean(disc_time.PACii(PACii).licks.pcorr(per_ii).group(grNo).data),'b','LineWidth', 3,'EdgeColor','none')
                    case 3
                        bar(bar_offset,mean(disc_time.PACii(PACii).licks.pcorr(per_ii).group(grNo).data),'y','LineWidth', 3,'EdgeColor','none')
                end
                
                try
                    CI = bootci(1000, {@mean, disc_time.PACii(PACii).licks.pcorr(per_ii).group(grNo).data},'type','cper');
                    plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                    plot(bar_offset*ones(1,length(disc_time.PACii(PACii).licks.pcorr(per_ii).group(grNo).data)),disc_time.PACii(PACii).licks.pcorr(per_ii).group(grNo).data,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
                catch
                end
                
            end
            bar_offset = bar_offset + 1;
            
        end
        
        ylim([0 4.5])
        
        title(['LDA decision time for licks'])
        
        
        xticks([1 2 3 5 6 7])
        xticklabels({'nwt', 'nH', 'nKO', 'pwt', 'pH', 'pKO'})
        
        
        ylabel('Decision time (sec)')
    end
    
    %Bar graph plot for peak discrimination time
    %Plot the average
    figNo = figNo +1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    
    set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
    hold on
    
    ax=gca;ax.LineWidth=3;
    
    bar_offset = 0;
    
    for per_ii=2:-1:1
        
        for grNo=1:max(handles_out.drgbchoices.group_no)
            bar_offset = bar_offset +1;
            
            
            switch grNo
                case 1
                    bar(bar_offset,mean(disc_time.PACii(PACii).peaks.pcorr(per_ii).group(grNo).data),'g','LineWidth', 3,'EdgeColor','none')
                case 2
                    bar(bar_offset,mean(disc_time.PACii(PACii).peaks.pcorr(per_ii).group(grNo).data),'b','LineWidth', 3,'EdgeColor','none')
                case 3
                    bar(bar_offset,mean(disc_time.PACii(PACii).peaks.pcorr(per_ii).group(grNo).data),'y','LineWidth', 3,'EdgeColor','none')
            end
            
            try
                CI = bootci(1000, {@mean, disc_time.PACii(PACii).peaks.pcorr(per_ii).group(grNo).data},'type','cper');
                plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                plot(bar_offset*ones(1,length(disc_time.PACii(PACii).peaks.pcorr(per_ii).group(grNo).data)),disc_time.PACii(PACii).peaks.pcorr(per_ii).group(grNo).data,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
            catch
            end
            
        end
        bar_offset = bar_offset + 1;
        
    end
    
    
    ylim([0 4.5])
    
    title(['LDA decision time for peak theta/' handles_out.drgbchoices.PACnames{PACii} ])
    
    
    xticks([1 2 3 5 6 7])
    xticklabels({'nwt', 'nH', 'nKO', 'pwt', 'pH', 'pKO'})
    
    
    ylabel('Decision time (sec)')
    
    %Plot for peak discrimination time vs lick discrimination time
    figNo = figNo +1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    
    set(hFig, 'units','normalized','position',[.2 .2 .4 .4])
    hold on
    
    ax=gca;ax.LineWidth=3;
    
    bar_offset = 0;
    
    per_ii=1;
    
    all_lick_times=[];
    all_peak_times=[];
    for grNo=1:max(handles_out.drgbchoices.group_no)
        
        
        switch grNo
            case 1
                plot(disc_time.PACii(PACii).licks.pcorr(per_ii).group(grNo).data,disc_time.PACii(PACii).peaks.pcorr(per_ii).group(grNo).data,'o','MarkerFaceColor','g','MarkerEdgeColor','g')
            case 2
                plot(disc_time.PACii(PACii).licks.pcorr(per_ii).group(grNo).data,disc_time.PACii(PACii).peaks.pcorr(per_ii).group(grNo).data,'o','MarkerFaceColor','b','MarkerEdgeColor','b')
            case 3
                plot(disc_time.PACii(PACii).licks.pcorr(per_ii).group(grNo).data,disc_time.PACii(PACii).peaks.pcorr(per_ii).group(grNo).data,'o','MarkerFaceColor','y','MarkerEdgeColor','y')
        end
        
         all_lick_times=[all_lick_times disc_time.PACii(PACii).licks.pcorr(per_ii).group(grNo).data];
         all_peak_times=[all_peak_times disc_time.PACii(PACii).peaks.pcorr(per_ii).group(grNo).data];
        
    end
    
    
    %     end
    
    plot([0 1],[0 1],'-k')
    xlim([0 1])
    ylim([0 1])
    
    title(['LDA decision time for peak vs licks theta/' handles_out.drgbchoices.PACnames{PACii} ])
    

    ylabel('Peak decision time (sec)')
    xlabel('Lick decision time (sec)')
    
    [h,p]=ttest(all_lick_times,all_peak_times)
    
    %Bar graph plot for trough discrimination time
    %Plot the average
    figNo = figNo +1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    
    set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
    hold on
    
    ax=gca;ax.LineWidth=3;
    
    bar_offset = 0;
    
    for per_ii=2:-1:1
        
        for grNo=1:max(handles_out.drgbchoices.group_no)
            bar_offset = bar_offset +1;
            
            
            switch grNo
                case 1
                    bar(bar_offset,mean(disc_time.PACii(PACii).trough.pcorr(per_ii).group(grNo).data),'g','LineWidth', 3,'EdgeColor','none')
                case 2
                    bar(bar_offset,mean(disc_time.PACii(PACii).trough.pcorr(per_ii).group(grNo).data),'b','LineWidth', 3,'EdgeColor','none')
                case 3
                    bar(bar_offset,mean(disc_time.PACii(PACii).trough.pcorr(per_ii).group(grNo).data),'y','LineWidth', 3,'EdgeColor','none')
            end
            
            try
                CI = bootci(1000, {@mean, disc_time.PACii(PACii).trough.pcorr(per_ii).group(grNo).data},'type','cper');
                plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                plot(bar_offset*ones(1,length(disc_time.PACii(PACii).trough.pcorr(per_ii).group(grNo).data)),disc_time.PACii(PACii).trough.pcorr(per_ii).group(grNo).data,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
            catch
            end
            
        end
        bar_offset = bar_offset + 1;
        
    end
    
    ylim([0 4.5])
    
    title(['LDA decision time for trough theta/' handles_out.drgbchoices.PACnames{PACii} ])
    
    
    xticks([1 2 3 5 6 7])
    xticklabels({'nwt', 'nH', 'nKO', 'pwt', 'pH', 'pKO'})
    
    
    ylabel('Decision time (sec)')
    
    %Plot for peak discrimination time vs lick discrimination time
    figNo = figNo +1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    
    set(hFig, 'units','normalized','position',[.2 .2 .4 .4])
    hold on
    
    ax=gca;ax.LineWidth=3;
    
    bar_offset = 0;
    
    per_ii=1;
    
    all_lick_times=[];
    all_trough_times=[];
    for grNo=1:max(handles_out.drgbchoices.group_no)
        
        
        switch grNo
            case 1
                plot(disc_time.PACii(PACii).licks.pcorr(per_ii).group(grNo).data,disc_time.PACii(PACii).trough.pcorr(per_ii).group(grNo).data,'o','MarkerFaceColor','g','MarkerEdgeColor','g')
            case 2
                plot(disc_time.PACii(PACii).licks.pcorr(per_ii).group(grNo).data,disc_time.PACii(PACii).trough.pcorr(per_ii).group(grNo).data,'o','MarkerFaceColor','b','MarkerEdgeColor','b')
            case 3
                plot(disc_time.PACii(PACii).licks.pcorr(per_ii).group(grNo).data,disc_time.PACii(PACii).trough.pcorr(per_ii).group(grNo).data,'o','MarkerFaceColor','y','MarkerEdgeColor','y')
        end
        
         all_lick_times=[all_lick_times disc_time.PACii(PACii).licks.pcorr(per_ii).group(grNo).data];
         all_trough_times=[all_trough_times disc_time.PACii(PACii).trough.pcorr(per_ii).group(grNo).data];
        
    end
    
    
    %     end
    
    plot([0 1],[0 1],'-k')
    xlim([0 1])
    ylim([0 1])
    
    title(['LDA decision time for trough vs licks theta/' handles_out.drgbchoices.PACnames{PACii} ])
     

    ylabel('Trough decision time (sec)')
    xlabel('Lick decision time (sec)')
    
    [h,p]=ttest(all_lick_times,all_trough_times)
    
    %Perform the glm
    fprintf(1, ['\n\nglm for lick decision time for Theta/' handles_out.drgbchoices.PACnames{PACii} '\n'])
    tbl = table(glm_d_time.data',glm_d_time.group',glm_d_time.perCorr',glm_d_time.ptl',...
        'VariableNames',{'mean_pc','group','perCorr','lick_trough_peak'});
    mdl = fitglm(tbl,'mean_pc~group+perCorr+lick_trough_peak+group*lick_trough_peak*perCorr'...
        ,'CategoricalVars',[2,3,4])
    
    %Do ranksum/t test
    fprintf(1, ['\n\nRanksum or t-test p values for area under the curve for Theta/' handles_out.drgbchoices.PACnames{PACii} ' ' handles_out.drgbchoices.group_no_names{groupNo}  ' ' handles_out.drgbchoices.per_lab{percent_correct_ii} '\n'])
    try
        [output_data] = drgMutiRanksumorTtest(t_det_stats);
        fprintf(1, '\n\n')
    catch
    end
    
    lick_out.PACii(PACii).disc_time.PACii(PACii)=disc_time.PACii(PACii);
    
    pffft=1;
end


%Finally plot the fraction of the time that the animal is licking

PACii=1;


for percent_correct_ii=1:2
    
    for groupNo=1:max(handles_out.drgbchoices.group_no)
        
        %Gather all the data
        no_mice=0;
        no_mice_included=0;
        all_discriminant_splus_lickf=zeros(max(handles_out.drgbchoices.mouse_no),length(t));
        all_discriminant_sminus_lickf=zeros(max(handles_out.drgbchoices.mouse_no),length(t));
        
        for mouseNo=1:length(handles_out.discriminant_PACwavepower)
            try
                if handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated==1
                    per_ii=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).no_trials;
                    no_mice=no_mice+1;
                    if (per_ii>=20)&(sum(no_mice==mice_excluded)==0)
                        no_mice_included=no_mice_included+1;
                        flick_out.PACii(PACii).pcorr(percent_correct_ii).group(groupNo).mouseNos(no_mice_included)=mouseNo;
                        sp_trials=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).which_events(2,:);
                        these_sp_lickf=zeros(sum(sp_trials),length(t));
                        these_sp_lickf(:,:)=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).all_licks_per_tPACwave(logical(sp_trials),:);
                        all_discriminant_splus_lickf(no_mice_included,:)=mean(these_sp_lickf,1);
                        sm_trials=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).which_events(5,:);
                        these_sm_lickf=zeros(sum(sm_trials),length(t));
                        these_sm_lickf(:,:)=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).all_licks_per_tPACwave(logical(sm_trials),:);
                        all_discriminant_sminus_lickf(no_mice_included,:)=mean(these_sm_lickf,1);
                    end
                end
            catch
            end
        end
        
         all_discriminant_splus_lickf=all_discriminant_splus_lickf(1:no_mice_included,:);
         all_discriminant_sminus_lickf=all_discriminant_sminus_lickf(1:no_mice_included,:);
         
        no_mice_per(percent_correct_ii)=no_mice_included;
        
        %Plot the fractional lick
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        hFig=figure(figNo);
        
        set(hFig, 'units','normalized','position',[.61 .4 .38 .38])
        
        hold on
        
        ax=gca;ax.LineWidth=3;
        
        %First plot Sminus
        mean_sminus_lickf=nanmean(all_discriminant_sminus_lickf,1)';
        if size(all_discriminant_sminus_lickf,1)>2
            CIlickf = bootci(1000, {@nanmean, all_discriminant_sminus_lickf})';
            CIlickf(:,1)=mean_sminus_lickf-CIlickf(:,1);
            CIlickf(:,2)=CIlickf(:,2)-mean_sminus_lickf;
            if percent_correct_ii==1
                [hlick, hpCR] = boundedline(t,mean_sminus_lickf, CIlickf, 'cmap',[0 114/255 178/255]);
            else
                [hlick, hpCR] = boundedline(t,mean_sminus_lickf, CIlickf, 'cmap',[80/255 194/255 255/255]);
            end
        else
            plot(t,mean_sminus_lickf,'-k')
        end
        
        data=[];
        for mouseNo=1:no_mice_included
            data=[data mean(all_discriminant_sminus_lickf(mouseNo,(t>=odor_from)&(t<=odor_to)),2)];
        end
        
        lickf_out.PACii(PACii).lick.pcorr(percent_correct_ii).group(groupNo).sminus_data=data;
        
        %Then plot Splus
        mean_splus_lickf=nanmean(all_discriminant_splus_lickf,1)';
        if size(all_discriminant_splus_lickf,1)>2
            CIlickf = bootci(1000, {@nanmean, all_discriminant_splus_lickf})';
            CIlickf(:,1)=mean_splus_lickf-CIlickf(:,1);
            CIlickf(:,2)=CIlickf(:,2)-mean_splus_lickf;
            if percent_correct_ii==1
                [hlick, hpCR] = boundedline(t,mean_splus_lickf, CIlickf, 'cmap',[158/255 31/255 99/255]);
            else
                [hlick, hpCR] = boundedline(t,mean_splus_lickf, CIlickf, 'cmap',[238/255 111/255 179/255]);
            end
        else
            plot(t,mean_splus_lickf,'-k')
        end
        
        data=[];
        for mouseNo=1:no_mice_included
            data=[data mean(all_discriminant_splus_lickf(mouseNo,(t>=odor_from)&(t<=odor_to)),2)];
        end
        
        lickf_out.PACii(PACii).lick.pcorr(percent_correct_ii).group(groupNo).splus_data=data;
        
        %Odor on markers
        
        odorhl=plot([0 2.5],[0 0],'-k','LineWidth',5);
        
        
        title(['Lick fraction  for ' handles_out.drgbchoices.group_no_names{groupNo}  ' ' handles_out.drgbchoices.per_lab{percent_correct_ii}])
        
        xlabel('Time (sec)')
        ylabel(['lick fraction'])
        
        
    end
end
 
    
output_name=[pname 'pcorr_' fname];
save(output_name,'pcorr_out','lick_out','lickf_out','pval_out','disc_time','-v7.3')
 
fprintf(1, ['Finished processing ' discriminant_name '\n'])

pffft=1;

           


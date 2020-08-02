function drgAnalyzeLFPDiscriminantCaMKII
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

PACnames{1}='Beta';
PACnames{2}='Low gamma';
PACnames{3}='High gamma';

GenoNames{1}='WT';
GenoNames{2}='Het';
GenoNames{3}='KO';

%Plot the summary wavelet power at peak vs power at trough


figNo=0;

for fileNo=1:handlesdrgb.drgbchoices.no_files
    
    for PACii=[1 3]
        p_correct_stats=[];
        ii_stats=0;
   
        glm_ii=0;
        glm_correct=[];
       
%         AUCpeakloc1=[];
%         AUCpeakloc2=[];
        
        pname=handlesdrgb.drgbchoices.PathName{fileNo};
        fname=handlesdrgb.drgbchoices.FileName{fileNo};
        
        discriminant_name=[pname fname];
        load(discriminant_name)
        
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        hFig=figure(figNo);
        
        set(hFig, 'units','normalized','position',[.1 .1 .5 .5])
        
        x_offset=1;
        
        %Plot average percent correct for the LDA for peak and trough for
        %wavelet power referenced to PAC phase
        t=handles_out.t_power;
        
        for groupNo=1:max(handles_out.drgbchoices.group_no)
            
            for percent_correct_ii=2:-1:1
                
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
                
                if (percent_correct_ii==1)&(PACii==1)
                    fprintf(1, ['The number of mice included in the LDA analysis for odor pair ' handlesdrgb.drgbchoices.odorpair{fileNo} ' is %d\n\n\n'], no_mice_included)
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
                
                %             ii_stats=ii_stats+1;
                %             p_correct_stats(ii_stats).data=data;
                %             p_correct_stats(ii_stats).description=[handlesdrgb.drgbchoices.odorpair{fileNo} ' ' handles_out.drgbchoices.per_lab{percent_correct_ii} ' shuffled peak location ' num2str(handlesdrgb.drgbchoices.locations(fileNo))];
                %             p_correct_stats(ii_stats).data_ii=1;
                %             p_correct_stats(ii_stats).PACii=PACii;
                %             p_correct_stats(ii_stats).per_corr_ii=percent_correct_ii;
                %             p_correct_stats(ii_stats).groupNo=groupNo;
                %             %
                %
                %             glm_correct.data(glm_ii+1:glm_ii+length(data))=data;
                %             glm_correct.PACii(glm_ii+1:glm_ii+length(data))=PACii;
                %             glm_correct.group(glm_ii+1:glm_ii+length(data))=groupNo;
                %             glm_correct.perCorr(glm_ii+1:glm_ii+length(data))=percent_correct_ii;
                %             glm_correct.shuffled(glm_ii+1:glm_ii+length(data))=1;
                %             glm_correct.peak(glm_ii+1:glm_ii+length(data))=1;
                %             glm_correct.locations(glm_ii+1:glm_ii+length(data))=handlesdrgb.drgbchoices.locations(fileNo);
                %             glm_ii=glm_ii+length(data);
                
                data=[];
                for mouseNo=1:no_mice_included
                    data=[data (mean(all_discriminant_correct_shuffled_trough(mouseNo,(t>=auc_from)&(t<=auc_to)),2)-50)/50];
                end
                
                %             ii_stats=ii_stats+1;
                %             p_correct_stats(ii_stats).data=data;
                %             p_correct_stats(ii_stats).description=[handlesdrgb.drgbchoices.odorpair{fileNo} ' ' handles_out.drgbchoices.per_lab{percent_correct_ii} ' shuffled trough location ' num2str(handlesdrgb.drgbchoices.locations(fileNo))];
                %             p_correct_stats(ii_stats).data_ii=2;
                %             p_correct_stats(ii_stats).PACii=PACii;
                %             p_correct_stats(ii_stats).per_corr_ii=percent_correct_ii;
                %             p_correct_stats(ii_stats).groupNo=groupNo;
                %
                %             glm_correct.data(glm_ii+1:glm_ii+length(data))=data;
                %             glm_correct.PACii(glm_ii+1:glm_ii+length(data))=PACii;
                %             glm_correct.group(glm_ii+1:glm_ii+length(data))=groupNo;
                %             glm_correct.perCorr(glm_ii+1:glm_ii+length(data))=percent_correct_ii;
                %             glm_correct.shuffled(glm_ii+1:glm_ii+length(data))=1;
                %             glm_correct.peak(glm_ii+1:glm_ii+length(data))=0;
                %             glm_correct.locations(glm_ii+1:glm_ii+length(data))=handlesdrgb.drgbchoices.locations(fileNo);
                %             glm_ii=glm_ii+length(data);
                
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
                AUCtrough=data;
                ii_stats=ii_stats+1;
                p_correct_stats(ii_stats).data=data;
                p_correct_stats(ii_stats).description=[handlesdrgb.drgbchoices.odorpair{fileNo} ' ' handles_out.drgbchoices.per_lab{percent_correct_ii} ' ' GenoNames{groupNo}];
                p_correct_stats(ii_stats).PACii=PACii;
                p_correct_stats(ii_stats).per_corr_ii=percent_correct_ii;
                p_correct_stats(ii_stats).groupNo=groupNo;
                %
                glm_correct.data(glm_ii+1:glm_ii+length(data))=data;
                glm_correct.PACii(glm_ii+1:glm_ii+length(data))=PACii;
                glm_correct.group(glm_ii+1:glm_ii+length(data))=groupNo;
                glm_correct.perCorr(glm_ii+1:glm_ii+length(data))=percent_correct_ii;
                glm_correct.shuffled(glm_ii+1:glm_ii+length(data))=0;
                glm_correct.peak(glm_ii+1:glm_ii+length(data))=0;
%                 glm_correct.locations(glm_ii+1:glm_ii+length(data))=handlesdrgb.drgbchoices.locations(fileNo);
                glm_ii=glm_ii+length(data);
                
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
                AUCpeak=data;
                ii_stats=ii_stats+1;
                p_correct_stats(ii_stats).data=data;
                p_correct_stats(ii_stats).description=[handlesdrgb.drgbchoices.odorpair{fileNo} ' ' handles_out.drgbchoices.per_lab{percent_correct_ii} ' ' GenoNames{groupNo}];
                p_correct_stats(ii_stats).data_ii=4;
                p_correct_stats(ii_stats).PACii=PACii;
                p_correct_stats(ii_stats).per_corr_ii=percent_correct_ii;
                p_correct_stats(ii_stats).groupNo=groupNo;
                %
                glm_correct.data(glm_ii+1:glm_ii+length(data))=data;
                glm_correct.PACii(glm_ii+1:glm_ii+length(data))=PACii;
                glm_correct.group(glm_ii+1:glm_ii+length(data))=groupNo;
                glm_correct.perCorr(glm_ii+1:glm_ii+length(data))=percent_correct_ii;
                glm_correct.shuffled(glm_ii+1:glm_ii+length(data))=0;
                glm_correct.peak(glm_ii+1:glm_ii+length(data))=1;
%                 glm_correct.locations(glm_ii+1:glm_ii+length(data))=handlesdrgb.drgbchoices.locations(fileNo);
                glm_ii=glm_ii+length(data);
                
%                 if percent_correct_ii==1
%                     if handlesdrgb.drgbchoices.locations(fileNo)==1
%                         AUCpeakloc1=[AUCpeakloc1 data];
%                     else
%                         AUCpeakloc2=[AUCpeakloc2 data];
%                     end
%                 end
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
                
                %Plot peak
                subplot(2,1,1)
                hold on
                if percent_correct_ii==1
                    bar(x_offset,mean(AUCpeak),'r')
                else
                    bar(x_offset,mean(AUCpeak),'b')
                end
                CI = bootci(1000, {@mean, AUCpeak},'type','cper');
                plot([x_offset x_offset],CI,'-k','LineWidth',3)
                plot(x_offset*ones(1,length(AUCpeak)),AUCpeak,'o','MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
                                
                %Plot trough
                subplot(2,1,2)
                hold on
                if percent_correct_ii==1
                    bar(x_offset,mean(AUCtrough),'r')
                else
                    bar(x_offset,mean(AUCtrough),'b')
                end
                CI = bootci(1000, {@mean, AUCtrough},'type','cper');
                plot([x_offset x_offset],CI,'-k','LineWidth',3)
                plot(x_offset*ones(1,length(AUCtrough)),AUCtrough,'o','MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
                                
                if percent_correct_ii==1
                    x_offset=x_offset+2;
                else
                    x_offset=x_offset+1;
                end
                pfft=1;
            end
            pfft=1;
        end
        suptitle(['Discriminant analysis performance for theta/' PACnames{PACii} ' ' handlesdrgb.drgbchoices.odorpair{fileNo}])
        
        %Perform the glm for percent correct
        fprintf(1, ['\n\nglm for percent correct decoding in the LDA for Theta/' PACnames{PACii} '\n'])
        tbl = table(glm_correct.data',glm_correct.perCorr',glm_correct.peak',glm_correct.group',...
            'VariableNames',{'LDApcorr','proficiency','peak_trough','genotype'});
        mdl = fitglm(tbl,'LDApcorr~proficiency+genotype+peak_trough+proficiency*peak_trough*genotype'...
            ,'CategoricalVars',[2,3 4])
        
        %Do ranksum/t test
        fprintf(1, ['\n\nRanksum or t-test p values for percent correct decoding in the LDA for Theta/' PACnames{PACii} '\n'])
        try
            [output_data] = drgMutiRanksumorTtest(p_correct_stats);
            fprintf(1, '\n\n')
        catch
        end
    end
    
    
    
end

pffft=1;

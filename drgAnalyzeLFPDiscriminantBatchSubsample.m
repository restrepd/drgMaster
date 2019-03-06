function drgAnalyzeLFPDiscriminantBatchSubsample
%Analyzes the linear discriminant analysis performed by drgLFPDiscriminantBatch
%Takes as a in input the 'Discriminant_*.mat' output file from drgLFPDiscriminantBatch
%Performs an analysis of the timecourse for percent correct for LDA and for
%the PCA

close all
clear all

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

%Plot the subsample percent correct for each mouse
figNo=0;
each_mean_per_corr=[];
each_mean_per_corr_sh=[];
no_each_mean_per_corr=zeros(max(handles_out.drgbchoices.group_no),length(handles_out.drgbchoices.lowF),length(handles_out.drgbchoices.per_lab));
for mouseNo=1:max(handles_out.drgbchoices.mouse_no)
    for groupNo=1:max(handles_out.drgbchoices.group_no)
        
        figNo=figNo+1
        try
            close(figNo)
        catch
        end
        
        hFig=figure(figNo);
        set(hFig, 'units','normalized','position',[.07 0.1 .85 .7])
        
        keep_figure=1;
        for bwii=1:length(handles_out.drgbchoices.lowF)
            
            for percent_correct_ii=1:length(handles_out.drgbchoices.per_lab)
                
                subplot(length(handles_out.drgbchoices.per_lab),length(handles_out.drgbchoices.lowF),bwii+length(handles_out.drgbchoices.lowF)*(percent_correct_ii-1))
                hold on
                
                if isfield(handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii),'bwii')
                    par_out_fine=1;
                    try
                        par_out=handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).ecomb_par_out;
                    catch
                        par_out_fine=0;
                    end
                    if par_out_fine==1
                        for elNo=1:par_out(1).no_elect
                            
                            mean_dcsh(elNo)=mean(par_out(elNo).discriminant_correct_shuffled(1:par_out(elNo).no_samples));
                            tempCIdcsh = bootci(1000, {@mean, par_out(elNo).discriminant_correct_shuffled(1:par_out(elNo).no_samples)})';
                            CIdcsh(elNo,1)=mean_dcsh(elNo)-tempCIdcsh(1);
                            CIdcsh(elNo,2)=tempCIdcsh(2)-mean_dcsh(elNo);
                            
                            
                            mean_dc(elNo)=mean(par_out(elNo).discriminant_correct(1:par_out(elNo).no_samples));
                            tempCIdc = bootci(1000, {@mean, par_out(elNo).discriminant_correct(1:par_out(elNo).no_samples)})';
                            CIdc(elNo,1)=mean_dc(elNo)-tempCIdc(1);
                            CIdc(elNo,2)=tempCIdc(2)-mean_dc(elNo);
                            
                        end
                        no_each_mean_per_corr(groupNo,bwii,percent_correct_ii)=no_each_mean_per_corr(groupNo,bwii,percent_correct_ii)+1;
                        each_mean_per_corr(groupNo,bwii,percent_correct_ii,no_each_mean_per_corr(groupNo,bwii,percent_correct_ii),1:16)=mean_dc;
                        each_mean_per_corr_sh(groupNo,bwii,percent_correct_ii,no_each_mean_per_corr(groupNo,bwii,percent_correct_ii),1:16)=mean_dcsh;
                        
                        [hlCR, hpCR] = boundedline(1:par_out(1).no_elect,mean_dcsh, CIdcsh, 'b');
                        [hlCR, hpCR] = boundedline(1:par_out(1).no_elect,mean_dc, CIdc, 'r');
                        
                        %Plot tetrodes here
                        elNo=4;
                        tetNo=0;
                        mean_pcorr=[];
                        for sampNo=1:par_out(elNo).no_timepoints:par_out(elNo).no_samples
                            if par_out(elNo).is_tetrode_per_sample(sampNo)==1
                                tetNo=tetNo+1;
                                mean_pcorr(tetNo)=mean(par_out(elNo).discriminant_correct(sampNo:sampNo+par_out(elNo).no_timepoints-1));
                            end
                        end
                        
                        
                        plot(elNo*ones(1,tetNo),mean_pcorr,'ok','MarkerSize',5)
                        
                        
                        %Plot single electrodes here
                        elNo=1;
                        sampNo=1;
                        mean_pcorr=[];
                        for electrode_no=1:16
                            mean_pcorr(electrode_no)=mean(par_out(elNo).discriminant_correct(sampNo:sampNo+par_out(elNo).no_timepoints-1));
                            sampNo=sampNo+par_out(elNo).no_timepoints;
                        end
                        
                        plot(elNo*ones(1,16),mean_pcorr,'ok','MarkerSize',5)
                        
                        
                        xlim([1 par_out(1).no_elect])
                        ylim([40 110])
                        
                        title([handles_out.drgbchoices.bwlabels{bwii} ])
                        
                        
                        ylabel(['% correct '  handles_out.drgbchoices.per_lab{percent_correct_ii}])
                        xlabel('No electrodes')
                    else
                        keep_figure=0;
                    end
                end
            end
        end
        suptitle(['LDA Mouse No ' num2str(mouseNo) ' Group: ' handles_out.drgbchoices.group_no_names{groupNo}])
        
        if keep_figure==0
            close(hFig)
            figNo=figNo-1;
        end
        pffft=1
    end
end

%Plot the average
 for groupNo=1:max(handles_out.drgbchoices.group_no)
        
        figNo=figNo+1
        try
            close(figNo)
        catch
        end
        
        hFig=figure(figNo);
        set(hFig, 'units','normalized','position',[.07 0.1 .85 .7])
        
        keep_figure=1;
        for bwii=1:length(handles_out.drgbchoices.lowF)
            
            for percent_correct_ii=1:length(handles_out.drgbchoices.per_lab)
                
                subplot(length(handles_out.drgbchoices.per_lab),length(handles_out.drgbchoices.lowF),bwii+length(handles_out.drgbchoices.lowF)*(percent_correct_ii-1))
                hold on
                
                if no_each_mean_per_corr(groupNo,bwii,percent_correct_ii)>1
                    
                    these_mean_per_corr=zeros(no_each_mean_per_corr(groupNo,bwii,percent_correct_ii),16);
                    these_mean_per_corr(:,:)=each_mean_per_corr(groupNo,bwii,percent_correct_ii,1:no_each_mean_per_corr(groupNo,bwii,percent_correct_ii),1:16);
                    mean_dc=mean(these_mean_per_corr,1);
                    tempCIdc = bootci(1000, {@mean, these_mean_per_corr})';
                    CIdc=zeros(16,2);
                    CIdc(:,1)=mean_dc'-tempCIdc(:,1);
                    CIdc(:,2)=tempCIdc(:,2)-mean_dc';
                    
                    these_mean_per_corr_sh=zeros(no_each_mean_per_corr(groupNo,bwii,percent_correct_ii),16);
                    these_mean_per_corr_sh(:,:)=each_mean_per_corr_sh(groupNo,bwii,percent_correct_ii,1:no_each_mean_per_corr(groupNo,bwii,percent_correct_ii),1:16);
                    mean_dc_sh=mean(these_mean_per_corr_sh,1);
                    tempCIdc_sh = bootci(1000, {@mean, these_mean_per_corr_sh})';
                    CIdc_sh=zeros(16,2);
                    CIdc_sh(:,1)=mean_dc_sh'-tempCIdc_sh(:,1);
                    CIdc_sh(:,2)=tempCIdc_sh(:,2)-mean_dc_sh';
                    
                    [hlCR, hpCR] = boundedline(1:16,mean_dcsh, CIdcsh, 'b');
                    [hlCR, hpCR] = boundedline(1:16,mean_dc, CIdc, 'r');
                    
                    title([handles_out.drgbchoices.bwlabels{bwii} ])
                        
                        
                        ylabel(['% correct '  handles_out.drgbchoices.per_lab{percent_correct_ii}])
                        xlabel('No electrodes')
                else
                    keep_figure=0;
                end
            end
        end
        suptitle(['Average for all mice, Group: ' handles_out.drgbchoices.group_no_names{groupNo}])
        
        if keep_figure==0
            close(hFig)
            figNo=figNo-1;
        end
 end
 

 
 %Now do the aurocs
 
 each_auROC=[];
 each_auROC_sh=[];
 no_auROC=zeros(max(handles_out.drgbchoices.group_no),length(handles_out.drgbchoices.lowF),length(handles_out.drgbchoices.per_lab));
 for mouseNo=1:max(handles_out.drgbchoices.mouse_no)
     for groupNo=1:max(handles_out.drgbchoices.group_no)
         
         figNo=figNo+1
         try
             close(figNo)
         catch
         end
         
         hFig=figure(figNo);
         set(hFig, 'units','normalized','position',[.07 0.1 .85 .7])
         
         keep_figure=1;
         for bwii=1:length(handles_out.drgbchoices.lowF)
             
             for percent_correct_ii=1:length(handles_out.drgbchoices.per_lab)
                 
                 subplot(length(handles_out.drgbchoices.per_lab),length(handles_out.drgbchoices.lowF),bwii+length(handles_out.drgbchoices.lowF)*(percent_correct_ii-1))
                 hold on
                 
                 if isfield(handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii),'bwii')
                     par_out_fine=1;
                     try
                         par_out=handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).ecomb_par_out;
                     catch
                         par_out_fine=0;
                     end
                     if par_out_fine==1
                         
                         
                         
                         for elNo=1:par_out(1).no_elect
                             mean_auROC(elNo)=mean(par_out(elNo).auROC(1:par_out(elNo).no_samples));
                             tempCIauROC = bootci(1000, {@mean, par_out(elNo).auROC(1:par_out(elNo).no_samples)})';
                             CIauROC(elNo,1)=mean_auROC(elNo)-tempCIauROC(1);
                             CIauROC(elNo,2)=tempCIauROC(2)-mean_auROC(elNo);
                         end
                         
                         no_auROC(groupNo,bwii,percent_correct_ii)=no_auROC(groupNo,bwii,percent_correct_ii)+1;
                         each_auROC(groupNo,bwii,percent_correct_ii,no_auROC(groupNo,bwii,percent_correct_ii),1:16)=mean_auROC;
                         each_auROC_sh(groupNo,bwii,percent_correct_ii,no_auROC(groupNo,bwii,percent_correct_ii),1:16)=mean_dcsh;
                         
                         
                         [hlCR, hpCR] = boundedline(1:par_out(1).no_elect,mean_auROC, CIauROC, 'r');
                         
                         xlim([1 par_out(1).no_elect])
                         ylim([0 0.5])
                         
                         xlabel('No electrodes')
                         ylabel('auROC')
                         
                         
                     else
                         keep_figure=0;
                     end
                 end
             end
         end
         suptitle(['LDA auROC Mouse No ' num2str(mouseNo) ' Group: ' handles_out.drgbchoices.group_no_names{groupNo}])
         
         if keep_figure==0
             close(hFig)
             figNo=figNo-1;
         end
         pffft=1
     end
 end

pffft=1;


function drgAnalyzeLFPDiscriminantMultiConc
%Analyzes the linear discriminant analysis performed by drgLFPDiscriminantBatch
%Takes as a in input the 'drgbDiscPar' file listing 'Discriminant_*.mat' output files
%from drgLFPDiscriminantBatch
%
%Performs summary analises for LDA for different concentrations

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

figNo=0;

which_display=3;
mice_excluded=[];

[choiceFileName,choiceBatchPathName] = uigetfile({'drgbDiscPar*.m'},'Select the .m file with all the choices for analysis');
fprintf(1, ['\ndrgAnalyzeLFPDiscriminantMultiConc run for ' choiceFileName '\n\n']);

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
auc_from=0.5;
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

perCorrLabels{1}='Proficient';
perCorrLabels{2}='Naive';

p_AUCloc_stats=[];
ii_AUCloc_stats=0;
glm_AUCloc=[];
glm_AUC_ii=0;

%Initialize the AUC arrays
max_mouse_no=-5;
for fileNo=1:handlesdrgb.drgbchoices.no_files
    pname=handlesdrgb.drgbchoices.PathName;
    fname=handlesdrgb.drgbchoices.FileName{fileNo};
    discriminant_name=[pname fname];
    load(discriminant_name)
    max_mouse_no=max([max_mouse_no max(handles_out.drgbchoices.mouse_no)]);
end

AUC_peak=zeros(handlesdrgb.drgbchoices.no_files,2,2,max_mouse_no);
AUC_trough=zeros(handlesdrgb.drgbchoices.no_files,2,2,max_mouse_no);
AUC_shuffled_peak=zeros(handlesdrgb.drgbchoices.no_files,2,2,max_mouse_no);
AUC_shuffled_trough=zeros(handlesdrgb.drgbchoices.no_files,2,2,max_mouse_no);
num_mice_included=zeros(handlesdrgb.drgbchoices.no_files,2,2);
delta_conc=zeros(1,handlesdrgb.drgbchoices.no_files);
events1=zeros(1,handlesdrgb.drgbchoices.no_files);
events2=zeros(1,handlesdrgb.drgbchoices.no_files);

if groupNo==1
    concs=[10 3.2 1 0.32 0.1 0.032];
else
    concs=[0.032 0.1 0.32 1 3.2 10];
end

for PACii=[1 3]
    p_correct_stats=[];
    p_correct_stats_within=[];
    p_correct_stats_delta=[];
    ii_stats=0;
    p_dim_stats=[];
    ii_dim_stats=0;
    glm_ii=0;
    glm_correct=[];
    glm_dim_ii=0;
    glm_dim=[];
    
    AUCpeakloc1=[];
    AUCpeakloc2=[];
    
    for fileNo=1:handlesdrgb.drgbchoices.no_files
        
        pname=handlesdrgb.drgbchoices.PathName;
        fname=handlesdrgb.drgbchoices.FileName{fileNo};
        
        discriminant_name=[pname fname];
        load(discriminant_name)
        
        %Plot average percent correct for the LDA for peak and trough for
        %wavelet power referenced to PAC phase
        t=handles_out.t_power;
        
        
        for percent_correct_ii=1:2
 
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
            
            %Suffled peak
            data=[];
            for mouseNo=1:no_mice_included
                data=[data (mean(all_discriminant_correct_shuffled_peak(mouseNo,(t>=auc_from)&(t<=auc_to)),2)-50)/50];
            end
            

            AUC_shuffled_peak(fileNo,PACii,percent_correct_ii,1:no_mice_included)=data;
            num_mice_included(fileNo,PACii,percent_correct_ii)=no_mice_included;
            delta_conc(fileNo,PACii)=handlesdrgb.drgbchoices.delta_conc(fileNo);
            events1(fileNo,PACii)=handlesdrgb.drgbchoices.event1(fileNo);
            events2(fileNo,PACii)=handlesdrgb.drgbchoices.event2(fileNo);
            
            ii_stats=ii_stats+1;
            p_correct_stats(ii_stats).data=data;
            if percent_correct_ii==1
                if handlesdrgb.drgbchoices.within(fileNo)==1
                    p_correct_stats(ii_stats).description=[handlesdrgb.drgbchoices.odorpair{fileNo} ' naive shuffled peak within, delta= ' num2str(handlesdrgb.drgbchoices.delta_conc(fileNo))];
                else
                    p_correct_stats(ii_stats).description=[handlesdrgb.drgbchoices.odorpair{fileNo} ' naive shuffled peak between, delta= ' num2str(handlesdrgb.drgbchoices.delta_conc(fileNo))];
                end
            else
                if handlesdrgb.drgbchoices.within(fileNo)==1
                    p_correct_stats(ii_stats).description=[handlesdrgb.drgbchoices.odorpair{fileNo} ' proficient shuffled peak within, delta= ' num2str(handlesdrgb.drgbchoices.delta_conc(fileNo))];
                else
                    p_correct_stats(ii_stats).description=[handlesdrgb.drgbchoices.odorpair{fileNo} ' proficient shuffled peak between, delta= ' num2str(handlesdrgb.drgbchoices.delta_conc(fileNo))];
                end
            end

            glm_correct.data(glm_ii+1:glm_ii+length(data))=data;
            glm_correct.within(glm_ii+1:glm_ii+length(data))=handlesdrgb.drgbchoices.within(fileNo);
            glm_correct.shuffled(glm_ii+1:glm_ii+length(data))=1;
            glm_correct.peak(glm_ii+1:glm_ii+length(data))=1;
            glm_correct.perCorr(glm_ii+1:glm_ii+length(data))=percent_correct_ii;
            glm_correct.delta(glm_ii+1:glm_ii+length(data))=handlesdrgb.drgbchoices.delta_conc(fileNo);
            glm_ii=glm_ii+length(data);
            
            
            %Suffled trough
            data=[];
            for mouseNo=1:no_mice_included
                data=[data (mean(all_discriminant_correct_shuffled_trough(mouseNo,(t>=auc_from)&(t<=auc_to)),2)-50)/50];
            end
            
            AUC_shuffled_trough(fileNo,PACii,percent_correct_ii,1:no_mice_included)=data;
            
            
            ii_stats=ii_stats+1;
            p_correct_stats(ii_stats).data=data;
            if percent_correct_ii==1
                if handlesdrgb.drgbchoices.within(fileNo)==1
                    p_correct_stats(ii_stats).description=[handlesdrgb.drgbchoices.odorpair{fileNo} ' naive shuffled trough within, delta= ' num2str(handlesdrgb.drgbchoices.delta_conc(fileNo))];
                else
                    p_correct_stats(ii_stats).description=[handlesdrgb.drgbchoices.odorpair{fileNo} ' naive shuffled torugh between, delta= ' num2str(handlesdrgb.drgbchoices.delta_conc(fileNo))];
                end
            else
                   if handlesdrgb.drgbchoices.within(fileNo)==1
                    p_correct_stats(ii_stats).description=[handlesdrgb.drgbchoices.odorpair{fileNo} ' proficient shuffled trough within, delta= ' num2str(handlesdrgb.drgbchoices.delta_conc(fileNo))];
                else
                    p_correct_stats(ii_stats).description=[handlesdrgb.drgbchoices.odorpair{fileNo} ' proficient shuffled torugh between, delta= ' num2str(handlesdrgb.drgbchoices.delta_conc(fileNo))];
                end
            end
            
            glm_correct.data(glm_ii+1:glm_ii+length(data))=data;
            glm_correct.within(glm_ii+1:glm_ii+length(data))=handlesdrgb.drgbchoices.within(fileNo);
            glm_correct.shuffled(glm_ii+1:glm_ii+length(data))=1;
            glm_correct.peak(glm_ii+1:glm_ii+length(data))=0;
            glm_correct.delta(glm_ii+1:glm_ii+length(data))=handlesdrgb.drgbchoices.delta_conc(fileNo);
            glm_correct.perCorr(glm_ii+1:glm_ii+length(data))=percent_correct_ii;
            glm_ii=glm_ii+length(data);
            
            
            %Now plot the AUC for the trough
            all_discriminant_correct_trough=all_discriminant_correct_trough(1:no_mice_included,:);
            mean_dc_trough=mean(all_discriminant_correct_trough,1)';
     
            data=[];
            for mouseNo=1:no_mice_included
                data=[data (mean(all_discriminant_correct_trough(mouseNo,(t>=auc_from)&(t<=auc_to)),2)-50)/50];
            end
            
            AUC_trough(fileNo,PACii,percent_correct_ii,1:no_mice_included)=data;
            
            ii_stats=ii_stats+1;
            p_correct_stats(ii_stats).data=data;
            
            
            if percent_correct_ii==1
                if handlesdrgb.drgbchoices.within(fileNo)==1
                    p_correct_stats(ii_stats).description=[handlesdrgb.drgbchoices.odorpair{fileNo} ' naive trhough within, delta= ' num2str(handlesdrgb.drgbchoices.delta_conc(fileNo))];
                else
                    p_correct_stats(ii_stats).description=[handlesdrgb.drgbchoices.odorpair{fileNo} ' naive trough between, delta= ' num2str(handlesdrgb.drgbchoices.delta_conc(fileNo))];
                end
            else
                if handlesdrgb.drgbchoices.within(fileNo)==1
                    p_correct_stats(ii_stats).description=[handlesdrgb.drgbchoices.odorpair{fileNo} ' proficient trhough within, delta= ' num2str(handlesdrgb.drgbchoices.delta_conc(fileNo))];
                else
                    p_correct_stats(ii_stats).description=[handlesdrgb.drgbchoices.odorpair{fileNo} ' proficient trough between, delta= ' num2str(handlesdrgb.drgbchoices.delta_conc(fileNo))];
                end
            end

            glm_correct.data(glm_ii+1:glm_ii+length(data))=data;
            glm_correct.within(glm_ii+1:glm_ii+length(data))=handlesdrgb.drgbchoices.within(fileNo);
            glm_correct.shuffled(glm_ii+1:glm_ii+length(data))=0;
            glm_correct.peak(glm_ii+1:glm_ii+length(data))=0;
            glm_correct.delta(glm_ii+1:glm_ii+length(data))=handlesdrgb.drgbchoices.delta_conc(fileNo);
            glm_correct.perCorr(glm_ii+1:glm_ii+length(data))=percent_correct_ii;
            glm_ii=glm_ii+length(data);
            
            %Now plot the AUC for the peak
            all_discriminant_correct_peak=all_discriminant_correct_peak(1:no_mice_included,:);
            mean_dc_peak=mean(all_discriminant_correct_peak,1)';
  
            data=[];
            for mouseNo=1:no_mice_included
                data=[data (mean(all_discriminant_correct_peak(mouseNo,(t>=auc_from)&(t<=auc_to)),2)-50)/50];
            end
         
            AUC_peak(fileNo,PACii,percent_correct_ii,1:no_mice_included)=data;
            
            ii_stats=ii_stats+1;
            p_correct_stats(ii_stats).data=data;
            
            if percent_correct_ii==1
                if handlesdrgb.drgbchoices.within(fileNo)==1
                    p_correct_stats(ii_stats).description=[handlesdrgb.drgbchoices.odorpair{fileNo} ' naive peak within, delta= ' num2str(handlesdrgb.drgbchoices.delta_conc(fileNo))];
                else
                    p_correct_stats(ii_stats).description=[handlesdrgb.drgbchoices.odorpair{fileNo} ' naive peak between, delta= ' num2str(handlesdrgb.drgbchoices.delta_conc(fileNo))];
                end
            else
                if handlesdrgb.drgbchoices.within(fileNo)==1
                    p_correct_stats(ii_stats).description=[handlesdrgb.drgbchoices.odorpair{fileNo} ' proficient peak within, delta= ' num2str(handlesdrgb.drgbchoices.delta_conc(fileNo))];
                else
                    p_correct_stats(ii_stats).description=[handlesdrgb.drgbchoices.odorpair{fileNo} ' proficient peak between, delta= ' num2str(handlesdrgb.drgbchoices.delta_conc(fileNo))];
                end
            end

            glm_correct.data(glm_ii+1:glm_ii+length(data))=data;
            glm_correct.within(glm_ii+1:glm_ii+length(data))=handlesdrgb.drgbchoices.within(fileNo);
            glm_correct.shuffled(glm_ii+1:glm_ii+length(data))=0;
            glm_correct.peak(glm_ii+1:glm_ii+length(data))=1;
            glm_correct.delta(glm_ii+1:glm_ii+length(data))=handlesdrgb.drgbchoices.delta_conc(fileNo);
            glm_correct.perCorr(glm_ii+1:glm_ii+length(data))=percent_correct_ii;
            glm_ii=glm_ii+length(data);
            

        
        end
        
        
        
    end
    
       %Perform the glm for percent correct
    fprintf(1, ['\n\nglm for decoding AUC for LDA for Theta/' PACnames{PACii} '\n'])
    tbl = table(glm_correct.data',glm_correct.perCorr',glm_correct.peak',glm_correct.delta',glm_correct.within',glm_correct.shuffled',...
        'VariableNames',{'AUC','proficiency','peak_trough','delta','within','shuffled'});
    mdl = fitglm(tbl,'AUC~proficiency+peak_trough+delta+within+shuffled+proficiency*peak_trough*delta*within*shuffled'...
        ,'CategoricalVars',[2,3 5,6])
    
    %Do ranksum/t test
    fprintf(1, ['\n\nRanksum or t-test p values for decoding AUC for the LDA for Theta/' PACnames{PACii} '\n'])
    try
        [output_data] = drgMutiRanksumorTtest(p_correct_stats);
        fprintf(1, '\n\n')
    catch
    end
    
    %Plot the AUC peak vs delta graph for peak
    for percent_correct_ii=1:2

        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        hFig=figure(figNo);
        set(hFig, 'units','normalized','position',[.3 .3 .3 .3])
        hold on
        
        for delta_ii=1:5
            for within_ii=0:1
                %Get all the AUCs
                thisAUCp=[];
                thisAUCp_sh=[];
                
                for fileNo=1:handlesdrgb.drgbchoices.no_files
                    if (handlesdrgb.drgbchoices.within(fileNo)==within_ii)&(handlesdrgb.drgbchoices.delta_conc(fileNo)==delta_ii)
                        num_m=num_mice_included(fileNo,PACii,percent_correct_ii);
                        
                        one_AUCp=zeros(1,num_m);
                        one_AUCp(1,:)=AUC_peak(fileNo,PACii,percent_correct_ii,1:num_m);
                        thisAUCp=[thisAUCp one_AUCp];
                        
                        one_AUCp_sh=zeros(1,num_m);
                        one_AUCp_sh(1,:)=AUC_shuffled_peak(fileNo,PACii,percent_correct_ii,1:num_m);
                        thisAUCp_sh=[thisAUCp_sh one_AUCp_sh];
                    end
                end
                
                if ~isempty(thisAUCp)
                    %Plot data
                    pffft=1;
                    mean_thisAUCp=mean(thisAUCp);
                    CI_thisAUCp = bootci(1000, {@mean, thisAUCp})';
                    
                    if within_ii==1
                        plot(delta_ii-0.1,mean_thisAUCp,'or','MarkerSize',12,'MarkerFaceColor','r')
                        plot([delta_ii-0.1 delta_ii-0.1],CI_thisAUCp,'-k')
                        plot((delta_ii-0.1)*ones(1,length(thisAUCp)),thisAUCp,'ok','MarkerFaceColor','r')
                    else
                        plot(delta_ii+0.1,mean_thisAUCp,'ob','MarkerSize',12,'MarkerFaceColor','b')
                        plot([delta_ii+0.1 delta_ii+0.1],CI_thisAUCp,'-k')
                        plot((delta_ii+0.1)*ones(1,length(thisAUCp)),thisAUCp,'ok','MarkerFaceColor','b')
                    end
                end
                    
                if ~isempty(thisAUCp_sh)
                    %Shuffled
                    mean_thisAUCp_sh=mean(thisAUCp_sh);
                    CI_thisAUCp_sh = bootci(1000, {@mean, thisAUCp_sh})';
                    
                    if within_ii==1
                        plot(delta_ii-0.2,mean_thisAUCp_sh,'o','MarkerSize',12,'MarkerFaceColor',[1 0.7 0.7],'MarkerEdgeColor','r')
                        plot([delta_ii-0.2 delta_ii-0.2],CI_thisAUCp_sh,'-k')
                        plot((delta_ii-0.2)*ones(1,length(thisAUCp_sh)),thisAUCp_sh,'ok','MarkerFaceColor',[1 0.7 0.7])
                    else
                        plot(delta_ii+0.2,mean_thisAUCp_sh,'o','MarkerSize',12,'MarkerFaceColor',[0.7  0.7 1],'MarkerEdgeColor','b')
                        plot([delta_ii+0.2 delta_ii+0.2],CI_thisAUCp_sh,'-k')
                        plot((delta_ii+0.2)*ones(1,length(thisAUCp_sh)),thisAUCp_sh,'ok','MarkerFaceColor',[0.7  0.7 1])
                    end 
                end
                
                pfft=1;
                
            end
        end
        title(['AUC peak for ' perCorrLabels{percent_correct_ii} ' for Theta/' PACnames{PACii}])
        xlabel('Difference in concentration steps')
        ylabel('AUC peak')
        ylim([-0.2 1])
    end
    
    %Plot the AUC trough vs delta graph for trough
    for percent_correct_ii=1:2

        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        hFig=figure(figNo);
        set(hFig, 'units','normalized','position',[.3 .3 .3 .3])
        hold on
        
        for delta_ii=1:5
            for within_ii=0:1
                %Get all the AUCs
                thisAUCp=[];
                thisAUCp_sh=[];
                
                for fileNo=1:handlesdrgb.drgbchoices.no_files
                    if (handlesdrgb.drgbchoices.within(fileNo)==within_ii)&(handlesdrgb.drgbchoices.delta_conc(fileNo)==delta_ii)
                        num_m=num_mice_included(fileNo,PACii,percent_correct_ii);
                        
                        one_AUCp=zeros(1,num_m);
                        one_AUCp(1,:)=AUC_trough(fileNo,PACii,percent_correct_ii,1:num_m);
                        thisAUCp=[thisAUCp one_AUCp];
                        
                        one_AUCp_sh=zeros(1,num_m);
                        one_AUCp_sh(1,:)=AUC_shuffled_trough(fileNo,PACii,percent_correct_ii,1:num_m);
                        thisAUCp_sh=[thisAUCp_sh one_AUCp_sh];
                    end
                end
                
                if ~isempty(thisAUCp)
                    %Plot data
                    pffft=1;
                    mean_thisAUCp=mean(thisAUCp);
                    CI_thisAUCp = bootci(1000, {@mean, thisAUCp})';
                    
                    if within_ii==1
                        plot(delta_ii-0.1,mean_thisAUCp,'or','MarkerSize',12,'MarkerFaceColor','r')
                        plot([delta_ii-0.1 delta_ii-0.1],CI_thisAUCp,'-k')
                        plot((delta_ii-0.1)*ones(1,length(thisAUCp)),thisAUCp,'ok','MarkerFaceColor','r')
                    else
                        plot(delta_ii+0.1,mean_thisAUCp,'ob','MarkerSize',12,'MarkerFaceColor','b')
                        plot([delta_ii+0.1 delta_ii+0.1],CI_thisAUCp,'-k')
                        plot((delta_ii+0.1)*ones(1,length(thisAUCp)),thisAUCp,'ok','MarkerFaceColor','b')
                    end
                end
                    
                if ~isempty(thisAUCp_sh)
                    %Shuffled
                    mean_thisAUCp_sh=mean(thisAUCp_sh);
                    CI_thisAUCp_sh = bootci(1000, {@mean, thisAUCp_sh})';
                    
                    if within_ii==1
                        plot(delta_ii-0.2,mean_thisAUCp_sh,'o','MarkerSize',12,'MarkerFaceColor',[1 0.7 0.7],'MarkerEdgeColor','r')
                        plot([delta_ii-0.2 delta_ii-0.2],CI_thisAUCp_sh,'-k')
                        plot((delta_ii-0.2)*ones(1,length(thisAUCp_sh)),thisAUCp_sh,'ok','MarkerFaceColor',[1 0.7 0.7])
                    else
                        plot(delta_ii+0.2,mean_thisAUCp_sh,'o','MarkerSize',12,'MarkerFaceColor',[0.7  0.7 1],'MarkerEdgeColor','b')
                        plot([delta_ii+0.2 delta_ii+0.2],CI_thisAUCp_sh,'-k')
                        plot((delta_ii+0.2)*ones(1,length(thisAUCp_sh)),thisAUCp_sh,'ok','MarkerFaceColor',[0.7  0.7 1])
                    end 
                end
                
                pfft=1;
                
            end
        end
        title(['AUC trough for ' perCorrLabels{percent_correct_ii} ' for Theta/' PACnames{PACii}])
        xlabel('Difference in concentration steps')
        ylabel('AUC trough')
        ylim([-0.2 1])
    end
    
     
end

pfft=1;

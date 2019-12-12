function drgAnalyzeLFPDiscriminantMultiBatchPerEvent
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

per_legend{1}='proficient';
per_legend{2}='naive';

%Plot the summary wavelet power at peak vs power at trough
figNo=0;



for PACii=[1 3]
    p_correct_stats=[];
    ii_stats=0;
    
    glm_ii=0;
    glm_correct=[];
    
    
    for fileNo=1:handlesdrgb.drgbchoices.no_files
        
        pname=handlesdrgb.drgbchoices.PathName{fileNo};
        fname=handlesdrgb.drgbchoices.FileName{fileNo};
        
        discriminant_name=[pname fname];
        load(discriminant_name)
        
        
        
        
        
        %Plot average percent correct for the LDA for peak and trough for
        %wavelet power referenced to PAC phase
        t=handles_out.t_power;
        
        if fileNo==1
            Hit_decoding_accuracy=zeros(handlesdrgb.drgbchoices.no_files,length(t));
            Miss_decoding_accuracy=zeros(handlesdrgb.drgbchoices.no_files,length(t));
            CR_decoding_accuracy=zeros(handlesdrgb.drgbchoices.no_files,length(t));
            FA_decoding_accuracy=zeros(handlesdrgb.drgbchoices.no_files,length(t));
            shuffled_decoding_accuracy=zeros(handlesdrgb.drgbchoices.no_files,length(t));
        end
        
        Hit_correct_peak=zeros(2,length(t));
        No_Hits_peak=zeros(2,length(t));
        
        Sp_correct_peak=zeros(2,length(t));
        No_Sp_peak=zeros(2,length(t));
        
        Miss_correct_peak=zeros(2,length(t));
        No_Miss_peak=zeros(2,length(t));
        
        CR_correct_peak=zeros(2,length(t));
        No_CR_peak=zeros(2,length(t));
        
        Sm_correct_peak=zeros(2,length(t));
        No_Sm_peak=zeros(2,length(t));
        
        FA_correct_peak=zeros(2,length(t));
        No_FA_peak=zeros(2,length(t));
        
        shuffled_correct_peak=zeros(2,length(t));
        No_shuffled_peak=zeros(2,length(t));
        
        
        for percent_correct_ii=1:1
            
            figure((PACii-1)*2+percent_correct_ii)
            
            %Gather all the data; I will only do peaks
            no_mice=0;
            no_mice_included=0;
            
            
            for mouseNo=1:length(handles_out.discriminant_PACwavepower)
                try
                    if handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated==1
                        per_ii=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).no_trials;
                        no_mice=no_mice+1;
                        if (per_ii>=20)&(sum(no_mice==mice_excluded)==0)
                            no_mice_included=no_mice_included+1;
                            which_events=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).which_events;
                            test_out=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).test_out_per_timepoint_peak;
                            sh_test_out=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).shuffled_out_per_timepoint_peak;
                            for trNo=1:per_ii
                                for evNo=1:6
                                    if which_events(evNo,trNo)==1
                                        for t_ii=1:length(t)
                                            switch evNo
                                                case 1 %Hit
                                                    Hit_correct_peak(percent_correct_ii,t_ii)=Hit_correct_peak(percent_correct_ii,t_ii)+test_out(1,trNo,t_ii);
                                                    shuffled_correct_peak(percent_correct_ii,t_ii)=shuffled_correct_peak(percent_correct_ii,t_ii)+sh_test_out(1,trNo,t_ii);
                                                    No_Hits_peak(percent_correct_ii,t_ii)=No_Hits_peak(percent_correct_ii,t_ii)+1;
                                                    No_shuffled_peak(percent_correct_ii,t_ii)=No_shuffled_peak(percent_correct_ii,t_ii)+1;
                                                case 2 %S+
                                                    Sp_correct_peak(percent_correct_ii,t_ii)=Sp_correct_peak(percent_correct_ii,t_ii)+test_out(1,trNo,t_ii);
                                                    %                                                     shuffled_Sp_correct_peak(percent_correct_ii,t_ii)=shuffled_Sp_correct_peak(percent_correct_ii,t_ii)+sh_test_out(1,trNo,t_ii);
                                                    No_Sp_peak(percent_correct_ii,t_ii)=No_Sp_peak(percent_correct_ii,t_ii)+1;
                                                case 3 %Miss
                                                    Miss_correct_peak(percent_correct_ii,t_ii)=Miss_correct_peak(percent_correct_ii,t_ii)+test_out(1,trNo,t_ii);
                                                    shuffled_correct_peak(percent_correct_ii,t_ii)=shuffled_correct_peak(percent_correct_ii,t_ii)+sh_test_out(1,trNo,t_ii);
                                                    No_Miss_peak(percent_correct_ii,t_ii)=No_Miss_peak(percent_correct_ii,t_ii)+1;
                                                    No_shuffled_peak(percent_correct_ii,t_ii)=No_shuffled_peak(percent_correct_ii,t_ii)+1;
                                                case 4 %CR
                                                    CR_correct_peak(percent_correct_ii,t_ii)=CR_correct_peak(percent_correct_ii,t_ii)+test_out(2,trNo,t_ii);
                                                    shuffled_correct_peak(percent_correct_ii,t_ii)=shuffled_correct_peak(percent_correct_ii,t_ii)+sh_test_out(2,trNo,t_ii);
                                                    No_CR_peak(percent_correct_ii,t_ii)=No_CR_peak(percent_correct_ii,t_ii)+1;
                                                    No_shuffled_peak(percent_correct_ii,t_ii)=No_shuffled_peak(percent_correct_ii,t_ii)+1;
                                                case 5 %S-
                                                    Sm_correct_peak(percent_correct_ii,t_ii)=Sm_correct_peak(percent_correct_ii,t_ii)+test_out(2,trNo,t_ii);
                                                    %                                                     shuffled_Sm_correct_peak(percent_correct_ii,t_ii)=shuffled_Sm_correct_peak(percent_correct_ii,t_ii)+sh_test_out(2,trNo,t_ii);
                                                    No_Sm_peak(percent_correct_ii,t_ii)=No_Sm_peak(percent_correct_ii,t_ii)+1;
                                                case 6 %FA
                                                    FA_correct_peak(percent_correct_ii,t_ii)=FA_correct_peak(percent_correct_ii,t_ii)+test_out(2,trNo,t_ii);
                                                    shuffled_correct_peak(percent_correct_ii,t_ii)=shuffled_correct_peak(percent_correct_ii,t_ii)+sh_test_out(2,trNo,t_ii);
                                                    No_FA_peak(percent_correct_ii,t_ii)=No_FA_peak(percent_correct_ii,t_ii)+1;
                                                    No_shuffled_peak(percent_correct_ii,t_ii)=No_shuffled_peak(percent_correct_ii,t_ii)+1;
                                            end
                                        end
                                    end
                                end
                            end
                            
                        end
                    end
                catch
                end
            end
            
            if (percent_correct_ii==1)&(PACii==1)
                fprintf(1, ['The number of mice included in the LDA analysis for odor pair ' handlesdrgb.drgbchoices.odorpair{fileNo} ' is %d\n\n\n'], no_mice_included)
            end
            
            no_mice_per(percent_correct_ii)=no_mice_included;
            
            %Plot the per event time course
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            hFig=figure(figNo);
            set(hFig, 'units','normalized','position',[.3 .3 .3 .3])
            hold on
            
            %Hits
            these_correct=zeros(1,length(t));
            these_correct(1,:)=Hit_correct_peak(percent_correct_ii,:);
            all_these=zeros(1,length(t));
            all_these(1,:)=No_Hits_peak(percent_correct_ii,:);
            plot(t,100*these_correct./all_these,'-r')
            
            Hit_decoding_accuracy(fileNo,:)=100*these_correct./all_these;
            
            glm_correct.data(glm_ii+1:glm_ii+length(t))=100*these_correct./all_these;
            glm_correct.time(glm_ii+1:glm_ii+length(t))=t;
            glm_correct.perCorr(glm_ii+1:glm_ii+length(t))=percent_correct_ii;
            glm_correct.event(glm_ii+1:glm_ii+length(t))=1;
            glm_correct.expt(glm_ii+1:glm_ii+length(t))=handlesdrgb.drgbchoices.location(fileNo);
            glm_ii=glm_ii+length(t);
            
            %Miss
            these_correct=zeros(1,length(t));
            these_correct(1,:)=Miss_correct_peak(percent_correct_ii,:);
            all_these=zeros(1,length(t));
            all_these(1,:)=No_Miss_peak(percent_correct_ii,:);
            plot(t,100*these_correct./all_these,'-c')
            
            Miss_decoding_accuracy(fileNo,:)=100*these_correct./all_these;
            
            glm_correct.data(glm_ii+1:glm_ii+length(t))=100*these_correct./all_these;
            glm_correct.time(glm_ii+1:glm_ii+length(t))=t;
            glm_correct.perCorr(glm_ii+1:glm_ii+length(t))=percent_correct_ii;
            glm_correct.event(glm_ii+1:glm_ii+length(t))=2;
            glm_correct.expt(glm_ii+1:glm_ii+length(t))=handlesdrgb.drgbchoices.location(fileNo);
            glm_ii=glm_ii+length(t);
            
            %CR
            these_correct=zeros(1,length(t));
            these_correct(1,:)=CR_correct_peak(percent_correct_ii,:);
            all_these=zeros(1,length(t));
            all_these(1,:)=No_CR_peak(percent_correct_ii,:);
            plot(t,100*these_correct./all_these,'-b')
            
            CR_decoding_accuracy(fileNo,:)=100*these_correct./all_these;
            
            glm_correct.data(glm_ii+1:glm_ii+length(t))=100*these_correct./all_these;
            glm_correct.time(glm_ii+1:glm_ii+length(t))=t;
            glm_correct.perCorr(glm_ii+1:glm_ii+length(t))=percent_correct_ii;
            glm_correct.event(glm_ii+1:glm_ii+length(t))=3;
            glm_correct.expt(glm_ii+1:glm_ii+length(t))=handlesdrgb.drgbchoices.location(fileNo);
            glm_ii=glm_ii+length(t);
            
            %FA
            these_correct=zeros(1,length(t));
            these_correct(1,:)=FA_correct_peak(percent_correct_ii,:);
            all_these=zeros(1,length(t));
            all_these(1,:)=No_FA_peak(percent_correct_ii,:);
            plot(t,100*these_correct./all_these,'-m')
            
            FA_decoding_accuracy(fileNo,:)=100*these_correct./all_these;
            
            glm_correct.data(glm_ii+1:glm_ii+length(t))=100*these_correct./all_these;
            glm_correct.time(glm_ii+1:glm_ii+length(t))=t;
            glm_correct.perCorr(glm_ii+1:glm_ii+length(t))=percent_correct_ii;
            glm_correct.event(glm_ii+1:glm_ii+length(t))=4;
            glm_correct.expt(glm_ii+1:glm_ii+length(t))=handlesdrgb.drgbchoices.location(fileNo);
            glm_ii=glm_ii+length(t);
            
            %Shuffled
            these_correct=zeros(1,length(t));
            these_correct(1,:)=shuffled_correct_peak(percent_correct_ii,:);
            all_these=zeros(1,length(t));
            all_these(1,:)=No_shuffled_peak(percent_correct_ii,:);
            plot(t,100*these_correct./all_these,'-k')
            
            shuffled_decoding_accuracy(fileNo,:)=100*these_correct./all_these;
            
            glm_correct.data(glm_ii+1:glm_ii+length(t))=100*these_correct./all_these;
            glm_correct.time(glm_ii+1:glm_ii+length(t))=t;
            glm_correct.perCorr(glm_ii+1:glm_ii+length(t))=percent_correct_ii;
            glm_correct.event(glm_ii+1:glm_ii+length(t))=5;
            glm_correct.expt(glm_ii+1:glm_ii+length(t))=handlesdrgb.drgbchoices.location(fileNo);
            glm_ii=glm_ii+length(t);
            
            title(['Decoding accuracy theta/' PACnames{PACii} ' ' handlesdrgb.drgbchoices.odorpair{fileNo} ' ' per_legend{percent_correct_ii}])
            ylabel('Decoding accuracy')
            xlabel('Time(sec)')
            ylim([0 110])
            
        end
        
    end
    
    %Plot the overall time course
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    set(hFig, 'units','normalized','position',[.3 .3 .3 .3])
    hold on
    
  
    
    %Miss
    mean_decoding=mean(Miss_decoding_accuracy,1)';
    CIda = bootci(1000, {@mean, Miss_decoding_accuracy})';
    CIda(:,1)=mean_decoding-CIda(:,1);
    CIda(:,2)=CIda(:,2)-mean_decoding;
    
    [hlCR, hpCR] = boundedline(t,mean_decoding, CIda, 'c');
    p2=plot(t,mean_decoding, 'c','LineWidth',2);
    
      %Shuffled
    mean_decoding=mean(shuffled_decoding_accuracy,1)';
    CIda = bootci(1000, {@mean, shuffled_decoding_accuracy})';
    CIda(:,1)=mean_decoding-CIda(:,1);
    CIda(:,2)=CIda(:,2)-mean_decoding;
    
    [hlCR, hpCR] = boundedline(t,mean_decoding, CIda, 'k');
    p1=plot(t,mean_decoding, 'k','LineWidth',2);
    
    %CR
    mean_decoding=mean(CR_decoding_accuracy,1)';
    CIda = bootci(1000, {@mean, CR_decoding_accuracy})';
    CIda(:,1)=mean_decoding-CIda(:,1);
    CIda(:,2)=CIda(:,2)-mean_decoding;
    
    [hlCR, hpCR] = boundedline(t,mean_decoding, CIda, 'b');
    p3=plot(t,mean_decoding, 'b','LineWidth',2);
    
    %FA
    mean_decoding=mean(FA_decoding_accuracy,1)';
    CIda = bootci(1000, {@mean, FA_decoding_accuracy})';
    CIda(:,1)=mean_decoding-CIda(:,1);
    CIda(:,2)=CIda(:,2)-mean_decoding;
    
    [hlCR, hpCR] = boundedline(t,mean_decoding, CIda, 'm');
    p4=plot(t,mean_decoding, 'm','LineWidth',2);
    
    %Hits
    mean_decoding=mean(Hit_decoding_accuracy,1)';
    CIda = bootci(1000, {@mean, Hit_decoding_accuracy})';
    CIda(:,1)=mean_decoding-CIda(:,1);
    CIda(:,2)=CIda(:,2)-mean_decoding;
    
    [hlCR, hpCR] = boundedline(t,mean_decoding, CIda, 'r');
    p5=plot(t,mean_decoding, 'r','LineWidth',2);
    
    %Odor on markers
    plot([0 2.5],[20 20],'-k', 'LineWidth', 4)
    
    title(['Decoding accuracy theta/' PACnames{PACii} ' '  per_legend{percent_correct_ii}])
    ylabel('Decoding accuracy')
    xlabel('Time(sec)')
    ylim([0 110])
    
    %     for expNo=1:2
%         %Plot the overall time course
%         figNo=figNo+1;
%         try
%             close(figNo)
%         catch
%         end
%         hFig=figure(figNo);
%         set(hFig, 'units','normalized','position',[.3 .3 .3 .3])
%         hold on
%         
%         %Shuffled
%         mean_decoding=mean(shuffled_decoding_accuracy(handlesdrgb.drgbchoices.location==expNo,:),1)';
%         CIda = bootci(1000, {@mean, shuffled_decoding_accuracy(handlesdrgb.drgbchoices.location==expNo,:)})';
%         CIda(:,1)=mean_decoding-CIda(:,1);
%         CIda(:,2)=CIda(:,2)-mean_decoding;
%         
%         [hlCR, hpCR] = boundedline(t,mean_decoding, CIda, 'k');
%         p1=plot(t,mean_decoding, 'k','LineWidth',2);
%         
%         %Miss
%         mean_decoding=mean(Miss_decoding_accuracy(handlesdrgb.drgbchoices.location==expNo,:),1)';
%         CIda = bootci(1000, {@mean, Miss_decoding_accuracy(handlesdrgb.drgbchoices.location==expNo,:)})';
%         CIda(:,1)=mean_decoding-CIda(:,1);
%         CIda(:,2)=CIda(:,2)-mean_decoding;
%         
%         [hlCR, hpCR] = boundedline(t,mean_decoding, CIda, 'c');
%         p1=plot(t,mean_decoding, 'c','LineWidth',2);
%         
%         %CR
%         mean_decoding=mean(CR_decoding_accuracy(handlesdrgb.drgbchoices.location==expNo,:),1)';
%         CIda = bootci(1000, {@mean, CR_decoding_accuracy(handlesdrgb.drgbchoices.location==expNo,:)})';
%         CIda(:,1)=mean_decoding-CIda(:,1);
%         CIda(:,2)=CIda(:,2)-mean_decoding;
%         
%         [hlCR, hpCR] = boundedline(t,mean_decoding, CIda, 'b');
%         p1=plot(t,mean_decoding, 'b','LineWidth',2);
%         
%         %FA
%         mean_decoding=mean(FA_decoding_accuracy(handlesdrgb.drgbchoices.location==expNo,:),1)';
%         CIda = bootci(1000, {@mean, FA_decoding_accuracy(handlesdrgb.drgbchoices.location==expNo,:)})';
%         CIda(:,1)=mean_decoding-CIda(:,1);
%         CIda(:,2)=CIda(:,2)-mean_decoding;
%         
%         [hlCR, hpCR] = boundedline(t,mean_decoding, CIda, 'r');
%         p1=plot(t,mean_decoding, 'r','LineWidth',2);
%         
%         %Hits
%         mean_decoding=mean(Hit_decoding_accuracy(handlesdrgb.drgbchoices.location==expNo,:),1)';
%         CIda = bootci(1000, {@mean, Hit_decoding_accuracy(handlesdrgb.drgbchoices.location==expNo,:)})';
%         CIda(:,1)=mean_decoding-CIda(:,1);
%         CIda(:,2)=CIda(:,2)-mean_decoding;
%         
%         [hlCR, hpCR] = boundedline(t,mean_decoding, CIda, 'r');
%         p1=plot(t,mean_decoding, 'r','LineWidth',2);
%     end
    
    %Perform the glm for percent correct
     fprintf(1, ['\n\nglm LDA accuracy vs events and experiments theta/' PACnames{PACii} '\n'])
    tbl = table(glm_correct.data',glm_correct.time',glm_correct.event',glm_correct.expt',...
        'VariableNames',{'LDAaccuracy','time','events','experiments'});
    mdl = fitglm(tbl,'LDAaccuracy~time+events+experiments+time*events*experiments'...
        ,'CategoricalVars',[3,4])
    
    pffft=1;
    
    
end

pffft=1;
            


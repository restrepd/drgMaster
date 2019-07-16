function drgAnalyzeDimsLFPDiscriminantAllMice
%Analyzes the linear discriminant analysis performed by drgLFPDiscriminantBatch
%Takes as a in input the 'Discriminant_*.mat' output file from drgLFPDiscriminantBatch
%Performs an analysis of the timecourse for percent correct for LDA and for
%the PCA



warning('off')
close all
clear all



handles_outp=[];

t_odor_arrival=0.1;

which_display=3;
mice_excluded=[];

[fname,pname,nCancel] = uigetfile({'Discriminant_*.mat'},'Select the all mouse discriminant output file ...');
if nCancel
    inputPath = [pname,fname];
    pnameStart = pname;
    %save('timeCorr_cfg.mat','pnameStart','-append');
else
    error('Cancelled')
end

discriminant_name_all=[pname fname];
load(discriminant_name_all)

handles_all=handles_out;


%Used to troubleshoot the IAMO JL
% mice_included=[1 3 4];

figNo=0;
%Plot average percent correct for the LDA for peak and trough for
%wavelet power referenced to PAC phase
t=handles_out.t_power;
groupNo=1;

for PACii=1:length(handles_out.drgbchoices.PACburstLowF)
    
        
        %Plot decision time reationship
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        hFig=figure(figNo);
        hold on
        
        p_disct_stats=[];
        ii_stats=0;
        
        glm_ii=0;
        glm_disct=[];
        
        for percent_correct_ii=2:-1:1
            
            
            dim_trough=handles_all.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).dimensionality_trough;
            dim_peak=handles_all.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).dimensionality_peak;
            
            
            subplot(1,2,3-percent_correct_ii)
            hold on
            plot(t,dim_trough,'-b')
            plot(t,dim_peak,'-r')
            
            ylabel('Time (sec)')
            if percent_correct_ii==1
                legend('Trough','Peak')
            end
            ylim([0 15])
            
        end
        
        suptitle(['Dimensionality all mouse for Theta/' handles_out.drgbchoices.PACnames{PACii} ])
        
       
%         
%         %Perform the glm
%         fprintf(1, ['\n\nglm for discrimination times for Theta/' handles_out.drgbchoices.PACnames{PACii} '\n'])
%         tbl = table(glm_disct.data',glm_disct.peak_trough',glm_disct.percent_correct',...
%             'VariableNames',{'t_detect','peak_trough','proficient_naive'});
%         mdl = fitglm(tbl,'t_detect~peak_trough+proficient_naive+peak_trough*proficient_naive'...
%             ,'CategoricalVars',[2,3])
%         
%         %Do ranksum/t test
%         fprintf(1, ['\n\nRanksum or t-test p values for discrimination times for Theta/' handles_out.drgbchoices.PACnames{PACii} '\n'])
%         try
%             [output_data] = drgMutiRanksumorTtest(p_disct_stats);
%             fprintf(1, '\n\n')
%         catch
%         end
%         
%      
    
end

pffft=1;
 



%drgCaMKII_disc_time_summary.m
%Compares discrimination times between licks and LDA for tPRP
close all
clear all

figNo=0;

handles_out.drgbchoices.PACnames{1}='Beta';
handles_out.drgbchoices.PACnames{2}='Gamma';

handles_out.drgbchoices.group_no_names{1}='WT';
handles_out.drgbchoices.group_no_names{2}='Het';
handles_out.drgbchoices.group_no_names{3}='KO';

path='/Users/restrepd/Documents/Projects/CaMKII_analysis/Discriminant new 80/';

disc_file{1}='pcorr_Discriminant_CaMKIIAPEB_disc_PRP80_hippo_04142021.mat';
disc_file_hipp_vs_pre(1)=1;
disc_file{2}='pcorr_Discriminant_CaMKIIAPEB_disc_PRP80_pre_04232021.mat';
disc_file_hipp_vs_pre(2)=2;
disc_file{3}='pcorr_Discriminant_CaMKIIEAPA_dis80c_hipp2_04212021.mat';
disc_file_hipp_vs_pre(3)=1;
disc_file{4}='pcorr_Discriminant_CaMKIIEAPA_dis80c_pre2_04212021.mat';
disc_file_hipp_vs_pre(4)=2;
disc_file{5}='pcorr_Discriminant_CaMKIIEBAP_disc_80_hipp04182021.mat';
disc_file_hipp_vs_pre(5)=1;
disc_file{6}='pcorr_Discriminant_CaMKIIEBAP_disc_80_pre04182021.mat';
disc_file_hipp_vs_pre(6)=2;
disc_file{7}='pcorr_Discriminant_CaMKIIPAEA_disc80_hipp2_04162021.mat';
disc_file_hipp_vs_pre(7)=1;
disc_file{8}='pcorr_Discriminant_CaMKIIPAEA_disc80_pre2_04162021.mat';
disc_file_hipp_vs_pre(8)=2;
disc_file{9}='pcorr_Discriminant_CaMKIIpz1eapa_disc_80_04142021_hipp.mat';
disc_file_hipp_vs_pre(9)=1;
disc_file{10}='pcorr_Discriminant_CaMKIIpz1eapa_disc_80_04142021_pre.mat';
disc_file_hipp_vs_pre(10)=2;
disc_file{11}='pcorr_Discriminant_CaMKIIpz1paea_disc_04112021_80hipp.mat';
disc_file_hipp_vs_pre(11)=1;
disc_file{12}='pcorr_Discriminant_CaMKIIpz1paea_disc_04112021_80pref.mat';
disc_file_hipp_vs_pre(12)=2;
disc_file{13}='pcorr_Discriminant_CaMKIIpzz1ethylace_disc_80hipp2_04082021.mat';
disc_file_hipp_vs_pre(13)=1;
disc_file{14}='pcorr_Discriminant_CaMKIIpzz1ethylace_disc_80pref2_04072021.mat';
disc_file_hipp_vs_pre(14)=2;
disc_file{15}='pcorr_Discriminant_CaMKIIpzz1paea_disc80_04052021_hipp.mat';
disc_file_hipp_vs_pre(15)=1;
disc_file{16}='pcorr_Discriminant_CaMKIIpzz1paea_disc80_04052021_pref.mat';
disc_file_hipp_vs_pre(16)=2;

if length(handles_out.drgbchoices.PACnames)==3
    these_PACii=[1 3];
else
    these_PACii=[1 2];
end


for ii_hipp_vs_pre=1:2
    
    
    for PACii=these_PACii
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
        
        for grNo=1:3
            all_lick_times(grNo).data=[];
            all_peak_times(grNo).data=[];
        end
        
        for grNo=1:3
            
            
            for fileNo=1:16
                if disc_file_hipp_vs_pre(fileNo)==ii_hipp_vs_pre
                    load([path disc_file{fileNo}])
                    for grNo=1:3
                        all_lick_times(grNo).data=[all_lick_times(grNo).data disc_time.PACii(PACii).licks.pcorr(per_ii).group(grNo).data];
                        all_peak_times(grNo).data=[all_peak_times(grNo).data disc_time.PACii(PACii).peaks.pcorr(per_ii).group(grNo).data];
                    end
                end
            end
        end
        
        edges=[0:0.1:1.2];
        rand_offset=0.8;
        
        bar_offset = 1;
        
        glm_dt_peak=[];
        glm_ii=0;
        
        id_ii=0;
        input_data=[];
        
        for grNo=1:3
            
            switch grNo
                case 1
                    %Plot decision for licks
                    bar(bar_offset,mean(all_lick_times(grNo).data),'FaceColor',[0.7 0.7 0.7],'LineWidth', 3,'EdgeColor','none')
                    %Violin plot
                    [mean_out, CIout]=drgViolinPoint(all_lick_times(grNo).data...
                        ,edges,bar_offset,rand_offset,'k','k',2);
                    
                    bar_offset = bar_offset +1;
                    
                    %Plot decision for tPRP
                    bar(bar_offset,mean(all_peak_times(grNo).data),'g','LineWidth', 3,'EdgeColor','none')
                    %Violin plot
                    [mean_out, CIout]=drgViolinPoint(all_peak_times(grNo).data...
                        ,edges,bar_offset,rand_offset,'k','k',2);
                case 2
                    %Plot decision for licks
                    bar(bar_offset,mean(all_lick_times(grNo).data),'FaceColor',[0.7 0.7 0.7],'LineWidth', 3,'EdgeColor','none')
                    %Violin plot
                    [mean_out, CIout]=drgViolinPoint(all_lick_times(grNo).data...
                        ,edges,bar_offset,rand_offset,'k','k',2);
                    
                    bar_offset = bar_offset +1;
                    
                    %Plot decision for tPRP
                    bar(bar_offset,mean(all_peak_times(grNo).data),'b','LineWidth', 3,'EdgeColor','none')
                    %Violin plot
                    [mean_out, CIout]=drgViolinPoint(all_peak_times(grNo).data...
                        ,edges,bar_offset,rand_offset,'k','k',2);
                case 3
                    %Plot decision for licks
                    bar(bar_offset,mean(all_lick_times(grNo).data),'FaceColor',[0.7 0.7 0.7],'LineWidth', 3,'EdgeColor','none')
                    %Violin plot
                    [mean_out, CIout]=drgViolinPoint(all_lick_times(grNo).data...
                        ,edges,bar_offset,rand_offset,'k','k',2);
                    
                    bar_offset = bar_offset +1;
                    
                    %Plot decision for tPRP
                    bar(bar_offset,mean(all_peak_times(grNo).data),'y','LineWidth', 3,'EdgeColor','none')
                    %Violin plot
                    [mean_out, CIout]=drgViolinPoint(all_peak_times(grNo).data...
                        ,edges,bar_offset,rand_offset,'k','k',2);
            end
            bar_offset = bar_offset +2;
             
            %Licks
            glm_dt_peak.data(glm_ii+1:glm_ii+length(all_lick_times(grNo).data))=all_lick_times(grNo).data;
            glm_dt_peak.lick_vs_peak(glm_ii+1:glm_ii+length(all_lick_times(grNo).data))=0;
            glm_dt_peak.group(glm_ii+1:glm_ii+length(all_lick_times(grNo).data))=grNo;
            glm_ii=glm_ii+length(all_lick_times(grNo).data);
            
            id_ii=id_ii+1;
            input_data(id_ii).data=all_lick_times(grNo).data;
            input_data(id_ii).description=[handles_out.drgbchoices.group_no_names{grNo} ' ' ...
                'licks'];
            
            %Peak
            glm_dt_peak.data(glm_ii+1:glm_ii+length(all_peak_times(grNo).data))=all_peak_times(grNo).data;
            glm_dt_peak.lick_vs_peak(glm_ii+1:glm_ii+length(all_peak_times(grNo).data))=0;
            glm_dt_peak.group(glm_ii+1:glm_ii+length(all_peak_times(grNo).data))=grNo;
            glm_ii=glm_ii+length(all_peak_times(grNo).data);
            
            id_ii=id_ii+1;
            input_data(id_ii).data=all_peak_times(grNo).data;
            input_data(id_ii).description=[handles_out.drgbchoices.group_no_names{grNo} ' ' ...
                'peak'];
            
        end
        
        
%         ylim([0 1])
        
        if ii_hipp_vs_pre==1
            title(['LDA decision time for peak vs licks for hippocampus for theta/' handles_out.drgbchoices.PACnames{PACii} ])
        else
            title(['LDA decision time for peak vs licks for prefrontal for theta/' handles_out.drgbchoices.PACnames{PACii} ])
        end
        
        ylabel('Decision time (sec)')

        %Perform the glm
        fprintf(1, ['glm for discrimination time for peak tPRP for PAC theta' handles_out.drgbchoices.PACnames{PACii} '\n'])
        tbl = table(glm_dt_peak.data',glm_dt_peak.group',glm_dt_peak.lick_vs_peak',...
            'VariableNames',{'detection_dt','group','lick_vs_peak'});
        mdl = fitglm(tbl,'detection_dt~group+lick_vs_peak+lick_vs_peak*group'...
            ,'CategoricalVars',[2,3])
        
        
        %Do the ranksum/t-test
        fprintf(1, ['\n\nRanksum or t-test p values for discrimination time for peak tPRP for PAC theta' handles_out.drgbchoices.PACnames{PACii} '\n'])
        [output_data] = drgMutiRanksumorTtest(input_data);
        
        fprintf(1, ['\n\n'])
        
        
        %Plot for trough discrimination time vs lick discrimination time
        
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
        
        for grNo=1:3
            all_lick_times(grNo).data=[];
            all_trough_times(grNo).data=[];
        end
        
         for grNo=1:3
            
            for fileNo=1:16
                if disc_file_hipp_vs_pre(fileNo)==ii_hipp_vs_pre
                    load([path disc_file{fileNo}])
                    for grNo=1:3
                        all_lick_times(grNo).data=[all_lick_times(grNo).data disc_time.PACii(PACii).licks.pcorr(per_ii).group(grNo).data];
                        all_trough_times(grNo).data=[all_trough_times(grNo).data disc_time.PACii(PACii).trough.pcorr(per_ii).group(grNo).data];
                    end
                end
            end
        end
        
        edges=[0:0.1:1.2];
        rand_offset=0.8;
        
        bar_offset = 1;
        
        glm_dt_trough=[];
        glm_ii=0;
        
        id_ii=0;
        input_data=[];
        
        for grNo=1:3
            
            switch grNo
                case 1
                    %Plot decision for licks
                    bar(bar_offset,mean(all_lick_times(grNo).data),'FaceColor',[0.7 0.7 0.7],'LineWidth', 3,'EdgeColor','none')
                    %Violin plot
                    [mean_out, CIout]=drgViolinPoint(all_lick_times(grNo).data...
                        ,edges,bar_offset,rand_offset,'k','k',2);
                    
                    bar_offset = bar_offset +1;
                    
                    %Plot decision for tPRP
                    bar(bar_offset,mean(all_trough_times(grNo).data),'g','LineWidth', 3,'EdgeColor','none')
                    %Violin plot
                    [mean_out, CIout]=drgViolinPoint(all_trough_times(grNo).data...
                        ,edges,bar_offset,rand_offset,'k','k',2);
                case 2
                    %Plot decision for licks
                    bar(bar_offset,mean(all_lick_times(grNo).data),'FaceColor',[0.7 0.7 0.7],'LineWidth', 3,'EdgeColor','none')
                    %Violin plot
                    [mean_out, CIout]=drgViolinPoint(all_lick_times(grNo).data...
                        ,edges,bar_offset,rand_offset,'k','k',2);
                    
                    bar_offset = bar_offset +1;
                    
                    %Plot decision for tPRP
                    bar(bar_offset,mean(all_trough_times(grNo).data),'b','LineWidth', 3,'EdgeColor','none')
                    %Violin plot
                    [mean_out, CIout]=drgViolinPoint(all_trough_times(grNo).data...
                        ,edges,bar_offset,rand_offset,'k','k',2);
                case 3
                    %Plot decision for licks
                    bar(bar_offset,mean(all_lick_times(grNo).data),'FaceColor',[0.7 0.7 0.7],'LineWidth', 3,'EdgeColor','none')
                    %Violin plot
                    [mean_out, CIout]=drgViolinPoint(all_lick_times(grNo).data...
                        ,edges,bar_offset,rand_offset,'k','k',2);
                    
                    bar_offset = bar_offset +1;
                    
                    %Plot decision for tPRP
                    bar(bar_offset,mean(all_trough_times(grNo).data),'y','LineWidth', 3,'EdgeColor','none')
                    %Violin plot
                    [mean_out, CIout]=drgViolinPoint(all_trough_times(grNo).data...
                        ,edges,bar_offset,rand_offset,'k','k',2);
            end
            bar_offset = bar_offset +2;
            
            %Licks
            glm_dt_trough.data(glm_ii+1:glm_ii+length(all_lick_times(grNo).data))=all_lick_times(grNo).data;
            glm_dt_trough.lick_vs_trough(glm_ii+1:glm_ii+length(all_lick_times(grNo).data))=0;
            glm_dt_trough.group(glm_ii+1:glm_ii+length(all_lick_times(grNo).data))=grNo;
            glm_ii=glm_ii+length(all_lick_times(grNo).data);
            
            id_ii=id_ii+1;
            input_data(id_ii).data=all_lick_times(grNo).data;
            input_data(id_ii).description=[handles_out.drgbchoices.group_no_names{grNo} ' ' ...
                'licks'];
            
            %Peak
            glm_dt_trough.data(glm_ii+1:glm_ii+length(all_trough_times(grNo).data))=all_trough_times(grNo).data;
            glm_dt_trough.lick_vs_trough(glm_ii+1:glm_ii+length(all_trough_times(grNo).data))=0;
            glm_dt_trough.group(glm_ii+1:glm_ii+length(all_trough_times(grNo).data))=grNo;
            glm_ii=glm_ii+length(all_trough_times(grNo).data);
            
            id_ii=id_ii+1;
            input_data(id_ii).data=all_trough_times(grNo).data;
            input_data(id_ii).description=[handles_out.drgbchoices.group_no_names{grNo} ' ' ...
                'trough'];
            
        end
        
        
%         ylim([0 1])
        
        if ii_hipp_vs_pre==1
            title(['LDA decision time for trough vs licks for hippocampus for theta/' handles_out.drgbchoices.PACnames{PACii} ])
        else
            title(['LDA decision time for trough vs licks for prefrontal for theta/' handles_out.drgbchoices.PACnames{PACii} ])
        end
        
        ylabel('Decision time (sec)')
      
 
        %Perform the glm
        if ii_hipp_vs_pre==1
            fprintf(1, ['glm for discrimination time for hippocampus trough tPRP for PAC theta' handles_out.drgbchoices.PACnames{PACii} '\n'])
        else
            fprintf(1, ['glm for discrimination time for prefrontal trough tPRP for PAC theta' handles_out.drgbchoices.PACnames{PACii} '\n'])
        end
        tbl = table(glm_dt_trough.data',glm_dt_trough.group',glm_dt_trough.lick_vs_trough',...
            'VariableNames',{'detection_dt','group','lick_vs_trough'});
        mdl = fitglm(tbl,'detection_dt~group+lick_vs_trough+lick_vs_trough*group'...
            ,'CategoricalVars',[2,3])
        
        
        %Do the ranksum/t-test
        if ii_hipp_vs_pre==1
            fprintf(1, ['\n\nRanksum or t-test p values for discrimination time for hippocampus trough tPRP for PAC theta' handles_out.drgbchoices.PACnames{PACii} '\n'])
        else
            fprintf(1, ['\n\nRanksum or t-test p values for discrimination time for prefrontal trough tPRP for PAC theta' handles_out.drgbchoices.PACnames{PACii} '\n'])
            
        end
        [output_data] = drgMutiRanksumorTtest(input_data);
        
        fprintf(1, ['\n\n'])
        
    end
    
end

function drgSummaryBatchImgCohCaMKIIWT
%Analyzes the linear discriminant analysis performed by drgLFPDiscriminantBatch
%Takes as a in input the 'drgbDiscPar' file listing 'Discriminant_*.mat' output files
%from drgLFPDiscriminantBatch
%
%Performs summary analises for LDA and  PAC analysis performed with case 19
%of drgAnalysisBatchLFP

warning('off')

close all hidden
clear all



%If you want statistics to be done with the value for each odor pair for
%each mouse make this variable 1
mouse_op=1;

bandwidth_names{1}='Theta';
bandwidth_names{2}='Beta';
bandwidth_names{3}='Low gamma';
bandwidth_names{4}='High gamma';

prof_naive_leg{1}='Proficient';
prof_naive_leg{2}='Naive';

group_legend{1}='WT';
group_legend{2}='Het';
group_legend{3}='KO';

evTypeLabels{1}='S+';
evTypeLabels{2}='S-';

peak_label{1}='Trough';
peak_label{2}='Peak';

%Location of files
% hippPathName='E:\CaMKIIpaper\datos sumarry\coherence\';
hippPathName='/Users/restrepd/Documents/Projects/CaMKII_analysis/Imaginary coherence/';
%
% %Files
% FileName{1}='CaMKIIacetocohe02012021_out.mat';
% FileName{2}='CaMKIIethylbenacetocohe2262021_out.mat';
% FileName{3}='CaMKIIPAEAcohe02082021_out.mat';
% FileName{4}='CaMKIIEAPAcohe2262021_out.mat';
% FileName{5}='CaMKIIpz1EAcohe02142021_out.mat';
% FileName{6}='CaMKIIPZ1PAEAcohe202102021_out.mat';
% FileName{7}='CaMKIIpzz1EAPAcohe02112021_out.mat';
% FileName{8}='CaMKIIpzz1propylacecohe02092021_out.mat';


FileName{1}='CaMKIIAPEBimgcohe01222022_out80.mat';
FileName{2}='spm_LFP_EBAP_imgcoh01282022_out80.mat';
FileName{3}='CaMKIIPAEAimagcohe01112022_out80';
FileName{4}='CaMKIIEAPAimagcohe01072022_out80.mat';
FileName{5}='CaMKIIpz1EAPAimagcohe01172022_out80.mat';
FileName{6}='CaMKIIpz1PAEAimagcohe01122022_out80.mat';
FileName{7}='spm_LFP_pzz1EAPA_imgcoh01252022_out80.mat';
FileName{8}='CaMKIIpzz1PAEAimgcoh01222021_out80.mat';

%Load the table of mouse numbers
%Note: This may need to be revised for PRP
mouse_no_table='/Users/restrepd/Documents/Projects/CaMKII_analysis/Reply_to_reviewers/camkii_mice_per_odor_pair_for_img_coh.xlsx';
T_mouse_no = readtable(mouse_no_table);


%Text file for statistical output
fileID = fopen([hippPathName 'drgSummaryBatchImgCohCaMKIIWT.txt'],'w');


%Load data
all_files=[];

for ii=1:length(FileName)
    load([hippPathName FileName{ii}])
    all_files(ii).handles_out=handles_out;
end

figNo=0;

%Now plot the average PRP for each electrode calculated per mouse
%(including all sessions for each mouse)
edges=[-0.5:0.05:0.5];
rand_offset=0.8;


for bwii=[1 2 4]    %for amplitude bandwidths (beta, low gamma, high gamma)
    
    glm_coh=[];
    glm_ii=0;
    
    id_ii=0;
    input_data=[];
    
    glm_coh_mm=[];
    glm_ii_mm=0;
    
    id_ii_mm=0;
    input_data_mm=[];
    
    %Plot the average
    figNo = figNo +1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    ax=gca;ax.LineWidth=3;
    
    set(hFig, 'units','normalized','position',[.1 .5 .3 .4])
    hold on
    
    ax=gca;ax.LineWidth=3;
    
    bar_lab_loc=[];
    no_ev_labels=0;
    ii_gr_included=0;
    bar_offset = 0;
    
    for evNo=1:2
        
        %         these_coh_per_Ev=[];
        %         ii_coh=0;
        coh_thisEv=[];
        mouse_no_thisEv=[];
        per_ii_thisEv=[];
        for per_ii=2:-1:1
            
            grNo=1;
            
            bar_offset = bar_offset + 1;
            
            %Get these coherence values
            these_coh=[];
            ii_coh=0;
         
            for ii=1:length(FileName)
                this_jj=[];
                these_mouse_no_per_op=T_mouse_no.mouse_no_per_op(T_mouse_no.odor_pair_no==ii);
                these_mouse_no=T_mouse_no.mouse_no(T_mouse_no.odor_pair_no==ii);
                for jj=1:all_files(ii).handles_out.dcoh_ii
                    
                    if all_files(ii).handles_out.dcoh_values(jj).pacii==bwii
                        if all_files(ii).handles_out.dcoh_values(jj).evNo==evNo
                            if all_files(ii).handles_out.dcoh_values(jj).per_ii==per_ii
                                if all_files(ii).handles_out.dcoh_values(jj).groupNo==grNo
                                    this_jj=jj;
                                end
                            end 
                        end
                    end
                end
                if ~isempty(this_jj)
                    
                        these_coh(ii_coh+1:ii_coh+length(all_files(ii).handles_out.dcoh_values(this_jj).dcoh_per_mouse))=all_files(ii).handles_out.dcoh_values(this_jj).dcoh_per_mouse;
                        ii_coh=ii_coh+length(all_files(ii).handles_out.dcoh_values(this_jj).dcoh_per_mouse);
                 
                    
                    %Assign to overall mouse numbers
                    these_dcoh_per_mouse=all_files(ii).handles_out.dcoh_values(this_jj).dcoh_per_mouse;
                    these_mouse_nos=all_files(ii).handles_out.dcoh_values(this_jj).mouseNo;
                    these_mouse_no_overall=[];
                    for ii_mouse=1:length(these_mouse_nos)
                        these_mouse_no_overall(ii_mouse)=these_mouse_no(find(these_mouse_nos(ii_mouse)==these_mouse_no_per_op));
                    end
                    
                    coh_thisEv=[coh_thisEv these_dcoh_per_mouse];
                    mouse_no_thisEv=[mouse_no_thisEv these_mouse_no_overall];
                    per_ii_thisEv=[per_ii_thisEv per_ii*ones(1,length(these_dcoh_per_mouse))];

                end
            end
            
              
            if evNo==2
                if per_ii==1
                    %S- Proficient
                    bar(bar_offset,mean(these_coh),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
%                     drgViolinPlot([bar_offset-1 bar_offset],[last_these_coh;these_coh],edges,rand_offset,'k','k',4,1);
                else
                    %S- Naive
                    bar(bar_offset,mean(these_coh),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
                    last_these_coh=these_coh;
                end
            else
                if per_ii==1
                    %S+ Proficient
                    bar(bar_offset,mean(these_coh),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
%                     drgViolinPlot([bar_offset-1 bar_offset],[last_these_coh;these_coh],edges,rand_offset,'k','k',4,1);
                else
                    %S+ naive
                    bar(bar_offset,mean(these_coh),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
                    last_these_coh=these_coh;
                end
            end
            
%             ii_these_coh=ii_these_coh+1;
%             these_coh_per_Ev(ii_these_coh,:)=these_coh;
            
            %Violin plot
            
%             [mean_out, CIout]=drgViolinPoint(these_coh,edges,bar_offset,rand_offset,'k','k',2);
            %             CI = bootci(1000, {@mean, these_coh},'type','cper');
            %             plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
            %             plot(bar_offset*ones(1,length(these_coh)),these_coh,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
            %
            
            
            glm_coh.data(glm_ii+1:glm_ii+length(these_coh))=these_coh;
            glm_coh.mouse_no(glm_ii+1:glm_ii+length(these_coh))=mouse_no_thisEv(per_ii_thisEv==per_ii);
            glm_coh.perCorr(glm_ii+1:glm_ii+length(these_coh))=per_ii*ones(1,length(these_coh));
            glm_coh.event(glm_ii+1:glm_ii+length(these_coh))=evNo*ones(1,length(these_coh));
            glm_ii=glm_ii+length(these_coh);
            
            id_ii=id_ii+1;
            input_data(id_ii).data=these_coh;
            input_data(id_ii).description=[ evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
             
        end
%         drgViolinPlot([bar_offset-1 bar_offset],these_coh_per_Ev,edges,rand_offset,'k','k',4,1);
        per_ii2_mean_coh_mm=[];
        per_ii1_mean_coh_mm=[];
        for mouseNo=1:max(mouse_no_thisEv)
            per_ii1_mean_coh_mm=[per_ii1_mean_coh_mm mean(coh_thisEv((mouse_no_thisEv==mouseNo)&(per_ii_thisEv==1)))];
            per_ii2_mean_coh_mm=[per_ii2_mean_coh_mm mean(coh_thisEv((mouse_no_thisEv==mouseNo)&(per_ii_thisEv==2)))];
        end
        
        %per_ii=1
        glm_coh_mm.data(glm_ii_mm+1:glm_ii_mm+length(per_ii1_mean_coh_mm))=per_ii1_mean_coh_mm;
        %                 glm_coh_mm.group(glm_ii_mm+1:glm_ii_mm+length(these_per_sig))=grNo*ones(1,length(per_ii1_mean_coh_mm));
        glm_coh_mm.perCorr(glm_ii_mm+1:glm_ii_mm+length(per_ii1_mean_coh_mm))=ones(1,length(per_ii1_mean_coh_mm));
        glm_coh_mm.event(glm_ii_mm+1:glm_ii_mm+length(per_ii1_mean_coh_mm))=evNo*ones(1,length(per_ii1_mean_coh_mm));
        glm_ii_mm=glm_ii_mm+length(per_ii1_mean_coh_mm);
        
        id_ii_mm=id_ii_mm+1;
        input_data_mm(id_ii_mm).data=per_ii1_mean_coh_mm;
        input_data_mm(id_ii_mm).description=[ evTypeLabels{evNo} ' ' prof_naive_leg{1}];
        
         %per_ii=1
        glm_coh_mm.data(glm_ii_mm+1:glm_ii_mm+length(per_ii2_mean_coh_mm))=per_ii2_mean_coh_mm;
        %                 glm_coh_mm.group(glm_ii_mm+1:glm_ii_mm+length(these_per_sig))=grNo*ones(1,length(per_ii2_mean_coh_mm));
        glm_coh_mm.perCorr(glm_ii_mm+1:glm_ii_mm+length(per_ii2_mean_coh_mm))=2*ones(1,length(per_ii2_mean_coh_mm));
        glm_coh_mm.event(glm_ii_mm+1:glm_ii_mm+length(per_ii2_mean_coh_mm))=evNo*ones(1,length(per_ii2_mean_coh_mm));
        glm_ii_mm=glm_ii_mm+length(per_ii2_mean_coh_mm);
        
        id_ii_mm=id_ii_mm+1;
        input_data_mm(id_ii_mm).data=per_ii2_mean_coh_mm;
        input_data_mm(id_ii_mm).description=[ evTypeLabels{evNo} ' ' prof_naive_leg{2}];
        
        for mouseNo=1:length(per_ii2_mean_coh_mm)
            plot([bar_offset-1 bar_offset],[per_ii2_mean_coh_mm(mouseNo) per_ii1_mean_coh_mm(mouseNo)],'-ok','MarkerSize',6,'LineWidth',2,'MarkerFaceColor',[0.6 0.6 0.6],'MarkerEdgeColor',[0.6 0.6 0.6],'Color',[0.6 0.6 0.6])
        end
        
        %Violin plot
        [mean_out, CIout]=drgViolinPoint(these_coh...
            ,edges,bar_offset,rand_offset,'k','k',3);
        
        %Violin plot
        [mean_out, CIout]=drgViolinPoint(last_these_coh...
            ,edges,bar_offset-1,rand_offset,'k','k',3);
        
        
        bar_offset = bar_offset + 1;
        
    end
    
    title(['Delta imaginary coherence for ' bandwidth_names{bwii}])
    ylim([-0.1 0.1])
    
    %     %Annotations identifying groups
    %     x_interval=0.8/ii_gr_included;
    %     for ii=1:ii_gr_included
    %         annotation('textbox',[0.7*x_interval+x_interval*(ii-1) 0.7 0.3 0.1],'String',handles_drgb.drgbchoices.group_no_names{ groups_included(ii)},'FitBoxToText','on');
    %     end
    %
    %     %Proficient/Naive annotations
    %     annotation('textbox',[0.15 0.8 0.3 0.1],'String','Proficient','FitBoxToText','on','Color','r','LineStyle','none');
    %     annotation('textbox',[0.15 0.75 0.3 0.1],'String','Naive','FitBoxToText','on','Color','b','LineStyle','none');
    
    plot([0 6],[0 0],'-k','LineWidth',2)
    xlim([0 6])
    
    xticks([1 2 4 5])
    xticklabels({'nS+', 'pS+','nS-', 'pS-'})
    
    ylabel('Delta imaginary coherence')
    
    
    %Perform the glm per mouse per odor pair
    fprintf(1, ['glm for delta imaginary coherence per mouse per odor pair for '  bandwidth_names{bwii} '\n'])
    fprintf(fileID, ['glm for delta imaginary coherence per mouse per odor pair for '  bandwidth_names{bwii} '\n']);
    
    
    tbl = table(glm_coh.data',glm_coh.perCorr',glm_coh.event',...
        'VariableNames',{'Img_coh','sp_vs_sm','naive_vs_proficient'});
    mdl = fitglm(tbl,'Img_coh~sp_vs_sm+naive_vs_proficient+sp_vs_sm*naive_vs_proficient'...
        ,'CategoricalVars',[2,3])
    
    txt = evalc('mdl');
    txt=regexp(txt,'<strong>','split');
    txt=cell2mat(txt);
    txt=regexp(txt,'</strong>','split');
    txt=cell2mat(txt);
    
    fprintf(fileID,'%s\n', txt);
    
    
    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for delta imaginary coherence per mouse per odor pair for ' bandwidth_names{bwii} '\n'])
    fprintf(fileID, ['\n\nRanksum or t-test p values for delta imaginary coherence per mouse per odor pair for ' bandwidth_names{bwii} '\n']);
    [output_data] = drgMutiRanksumorTtest(input_data, fileID);
    
    %Nested ANOVAN
    %https://www.mathworks.com/matlabcentral/answers/491365-within-between-subjects-in-anovan
    nesting=[0 0 0; ... % This line indicates that group factor is not nested in any other factor.
         0 0 0; ... % This line indicates that perCorr is not nested in any other factor.
         1 1 0];    % This line indicates that mouse_no (the last factor) is nested under group, perCorr and event
                    % (the 1 in position 1 on the line indicates nesting under the first factor).
    figNo=figNo+1;
                     
    [p anovanTbl stats]=anovan(glm_coh.data,{glm_coh.perCorr glm_coh.event glm_coh.mouse_no},...
    'model','interaction',...
    'nested',nesting,...
    'varnames',{'naive_vs_proficient', 'sp_vs_sm','mouse_no'});
  
    fprintf(fileID, ['\n\nNested ANOVAN for delta imaginary coherence per mouse per odor pair for '  bandwidth_names{bwii} '\n'])
    drgWriteANOVANtbl(anovanTbl,fileID);
    fprintf(fileID, '\n\n');
    
    
    %Perform the glm per mouse
    fprintf(1, ['glm for delta imaginary coherence per mouse for '  bandwidth_names{bwii} '\n'])
%     fprintf(fileID, ['glm for delta imaginary coherence per mouse for '  bandwidth_names{bwii} '\n']);
    
    
    tbl = table(glm_coh_mm.data',glm_coh_mm.perCorr',glm_coh_mm.event',...
        'VariableNames',{'Img_coh','sp_vs_sm','naive_vs_proficient'});
    mdl = fitglm(tbl,'Img_coh~sp_vs_sm+naive_vs_proficient+sp_vs_sm*naive_vs_proficient'...
        ,'CategoricalVars',[2,3])
    
%     txt = evalc('mdl');
%     txt=regexp(txt,'<strong>','split');
%     txt=cell2mat(txt);
%     txt=regexp(txt,'</strong>','split');
%     txt=cell2mat(txt);
%     
%     fprintf(fileID,'%s\n', txt);
    
    
    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for delta imaginary coherence per mouse for ' bandwidth_names{bwii} '\n'])
%     fprintf(fileID, ['\n\nRanksum or t-test p values for delta imaginary coherence per mouse for ' bandwidth_names{bwii} '\n']);
    [output_data] = drgMutiRanksumorTtest(input_data_mm);
    
    
end

    

%Plot the bounded lines
maxlP=-200000;
minlP=200000;

frequency=all_files(1).handles_out.frequency;

figNo = figNo + 1;
try
    close(figNo)
catch
end
hFig=figure(figNo);

set(hFig, 'units','normalized','position',[.1 .5 .3 .4])


set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 3)
hold on

for evNo=1:2
    
    
    for per_ii=2:-1:1      %performance bins. blue = naive, red = proficient
        

        
        
        grNo=1;
        
        %Get the deltaCxy_af
        %Get these MI values
        these_deltaCxy_af=[];
        ii_coh=0;
        mice_per_odor_pair=zeros(1,length(FileName));
        for ii=1:length(FileName)
            this_jj=[];
            for jj=1:all_files(ii).handles_out.dcohaf_ii
                if all_files(ii).handles_out.dcohaf_values(jj).evNo==evNo
                    if all_files(ii).handles_out.dcohaf_values(jj).per_ii==per_ii
                        if all_files(ii).handles_out.dcohaf_values(jj).groupNo==grNo
                            this_jj=jj;
                        end
                    end
                end
            end
            if (~isempty(this_jj))&(~isempty(all_files(ii).handles_out.dcohaf_values(this_jj).dcoh_per_mouse)) 
                if mouse_op==1
                    sz_dcoh=size(all_files(ii).handles_out.dcohaf_values(this_jj).dcoh_per_mouse);
                    these_deltaCxy_af(ii_coh+1:ii_coh+sz_dcoh(1),1:sz_dcoh(2))=all_files(ii).handles_out.dcohaf_values(this_jj).dcoh_per_mouse;
                    mice_per_odor_pair(ii)=max([mice_per_odor_pair(ii),length(all_files(ii).handles_out.dcohaf_values(this_jj).mouseNo)]);
                    
                    ii_coh=ii_coh+sz_dcoh(1);
                else
                    ii_coh=ii_coh+1;
                    these_deltaCxy_af(ii_coh,:)=all_files(ii).handles_out.dcohaf_values(this_jj).dcohaf;
                end
            end
            fprintf(1, ['\nNumber of mice for file No %d is %d\n'],ii,length(all_files(ii).handles_out.dcohaf_values(this_jj).mouseNo))
        end
        
        
        
         
        
        mean_deltaCxy=[];
        mean_deltaCxy=mean(these_deltaCxy_af,1);
        
        CI=[];
        CI = bootci(1000, {@mean, these_deltaCxy_af})';
        maxlP=max([maxlP max(CI(:))]);
        minlP=min([minlP min(CI(:))]);
        CI(:,1)= mean_deltaCxy'-CI(:,1);
        CI(:,2)=CI(:,2)- mean_deltaCxy';
        
        
        
        
        if evNo==2
            if per_ii==1
                %S- Proficient
                [hlCR, hpCR] = boundedline(frequency',mean_deltaCxy', CI, 'cmap',[158/255 31/255 99/255]);
                %                     bar(bar_offset,mean(these_coh),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
            else
                %S- Naive
                [hlCR, hpCR] = boundedline(frequency',mean_deltaCxy', CI, 'cmap',[238/255 111/255 179/255]);
                %                     bar(bar_offset,mean(these_coh),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
            end
        else
            if per_ii==1
                %S+ Proficient
                [hlCR, hpCR] = boundedline(frequency',mean_deltaCxy', CI, 'cmap',[0 114/255 178/255]);
                %                     bar(bar_offset,mean(these_coh),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
            else
                %S+ naive
                [hlCR, hpCR] = boundedline(frequency',mean_deltaCxy', CI, 'cmap',[80/255 194/255 255/255]);
                %                     bar(bar_offset,mean(these_coh),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
            end
        end
        
        pffft=1;
        
      
    end
    
end

title(['delta coherence '])
xlabel('Frequency (Hz')
ylabel('delta coherence')


%Now plot the average Cxy during the odor epoch
edges=[-0.5:0.05:0.5];
rand_offset=0.8;


for bwii=[1 2 4]    %for amplitude bandwidths (beta, low gamma, high gamma)
    
    glm_coh=[];
    glm_ii=0;
    
    id_ii=0;
    input_data=[];
    
    %Plot the average
    figNo = figNo +1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    ax=gca;ax.LineWidth=3;
    
    set(hFig, 'units','normalized','position',[.1 .5 .3 .4])
    hold on
    
    ax=gca;ax.LineWidth=3;
    
    bar_lab_loc=[];
    no_ev_labels=0;
    ii_gr_included=0;
    bar_offset = 0;
    
    for evNo=1:2
        
        these_coh_per_Ev=[];
        ii_coh=0;
        
        for per_ii=2:-1:1
            
            grNo=1;
            
            bar_offset = bar_offset + 1;
            
            %Get these coherence values
            these_coh=[];
            ii_these_coh=0;
            for ii=1:length(FileName)
                this_jj=[];
                for jj=1:all_files(ii).handles_out.dcoh_ii
                    
                    if all_files(ii).handles_out.odor_Cxy_values(jj).bwii==bwii
                        if all_files(ii).handles_out.odor_Cxy_values(jj).evNo==evNo
                            if all_files(ii).handles_out.odor_Cxy_values(jj).per_ii==per_ii
                                
                                
                                if all_files(ii).handles_out.odor_Cxy_values(jj).grNo==grNo
                                    this_jj=jj;
                                end
                                
                            end
                            
                        end
                    end
                end
                if ~isempty(this_jj)
                    if mouse_op==1
                        these_coh(ii_coh+1:ii_coh+length(all_files(ii).handles_out.odor_Cxy_values(this_jj).mean_Cxy_per_mouse))=all_files(ii).handles_out.odor_Cxy_values(this_jj).mean_Cxy_per_mouse;
                        ii_coh=ii_coh+length(all_files(ii).handles_out.odor_Cxy_values(this_jj).mean_Cxy_per_mouse);
                    else
                        ii_coh=ii_coh+1;
                        these_coh(ii_coh)=all_files(ii).handles_out.odor_Cxy_values(this_jj).dcoh;
                    end
                    
                end
            end
            
              
            if evNo==2
                if per_ii==1
                    %S- Proficient
                    bar(bar_offset,mean(these_coh),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
                else
                    %S- Naive
                    bar(bar_offset,mean(these_coh),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
                end
            else
                if per_ii==1
                    %S+ Proficient
                    bar(bar_offset,mean(these_coh),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
                else
                    %S+ naive
                    bar(bar_offset,mean(these_coh),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
                end
            end
            
%             ii_these_coh=ii_these_coh+1;
%             these_coh_per_Ev(ii_these_coh,:)=these_coh;
            
            %Violin plot
            
            [mean_out, CIout]=drgViolinPoint(these_coh,edges,bar_offset,rand_offset,'k','k',3);
            %             CI = bootci(1000, {@mean, these_coh},'type','cper');
            %             plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
            %             plot(bar_offset*ones(1,length(these_coh)),these_coh,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
            %
            
            
            glm_coh.data(glm_ii+1:glm_ii+length(these_coh))=these_coh;
            %                 glm_coh.group(glm_ii+1:glm_ii+length(these_coh))=grNo*ones(1,length(these_coh));
            glm_coh.perCorr(glm_ii+1:glm_ii+length(these_coh))=per_ii*ones(1,length(these_coh));
            glm_coh.event(glm_ii+1:glm_ii+length(these_coh))=evNo*ones(1,length(these_coh));
            glm_ii=glm_ii+length(these_coh);
            
            id_ii=id_ii+1;
            input_data(id_ii).data=these_coh;
            input_data(id_ii).description=[ evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
             
        end
%         drgViolinPlot([bar_offset-1 bar_offset],these_coh_per_Ev,edges,rand_offset,'k','k',4,1);
        
        bar_offset = bar_offset + 1;
        
    end
    
    title(['Odor coherence for ' bandwidth_names{bwii}])
    ylim([-0.1 0.1])
    
    %     %Annotations identifying groups
    %     x_interval=0.8/ii_gr_included;
    %     for ii=1:ii_gr_included
    %         annotation('textbox',[0.7*x_interval+x_interval*(ii-1) 0.7 0.3 0.1],'String',handles_drgb.drgbchoices.group_no_names{ groups_included(ii)},'FitBoxToText','on');
    %     end
    %
    %     %Proficient/Naive annotations
    %     annotation('textbox',[0.15 0.8 0.3 0.1],'String','Proficient','FitBoxToText','on','Color','r','LineStyle','none');
    %     annotation('textbox',[0.15 0.75 0.3 0.1],'String','Naive','FitBoxToText','on','Color','b','LineStyle','none');
    
    plot([0 6],[0 0],'-k','LineWidth',2)
    xlim([0 6])
    
    xticks([1 2 4 5])
    xticklabels({'nS+', 'pS+','nS-', 'pS-'})
    
    ylabel('Delta coherence')
    
    
    %Perform the glm
    fprintf(1, ['glm for odor coherence per mouse per odor pair for '  bandwidth_names{bwii} '\n'])
%     fprintf(fileID, ['glm for odor coherence per mouse per odor pair for '  bandwidth_names{bwii} '\n']);
    
    tbl = table(glm_coh.data',glm_coh.perCorr',glm_coh.event',...
        'VariableNames',{'odor_coh','perCorr','event'});
    mdl = fitglm(tbl,'odor_coh~perCorr+event+perCorr*event'...
        ,'CategoricalVars',[2,3]);
    
    txt = evalc('mdl');
    txt=regexp(txt,'<strong>','split');
    txt=cell2mat(txt);
    txt=regexp(txt,'</strong>','split');
    txt=cell2mat(txt);
    
%     fprintf(fileID,'%s\n', txt);
    
    
    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for odor coherence per mouse per odor pair for ' bandwidth_names{bwii} ' \n'])
%     fprintf(fileID, ['\n\nRanksum or t-test p values for odor coherence per mouse per odor pair for ' bandwidth_names{bwii} ' \n']);
    [output_data] = drgMutiRanksumorTtest(input_data);
    
    
end

%Now plot the cumulative histos
edges=[-0.5:0.05:0.5];
rand_offset=0.8;


for bwii=[1 2 4]    %for amplitude bandwidths (beta, low gamma, high gamma)
    
    glm_coh=[];
    glm_ii=0;
    
    id_ii=0;
    input_data=[];
    

    
    ax=gca;ax.LineWidth=3;
    
    bar_lab_loc=[];
    no_ev_labels=0;
    ii_gr_included=0;
    bar_offset = 0;
    
    for evNo=1:2
        
            %Plot the average
    figNo = figNo +1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    ax=gca;ax.LineWidth=3;
    
    set(hFig, 'units','normalized','position',[.1 .5 .3 .4])
    hold on
    
        these_coh_per_Ev=[];
        ii_coh=0;
        
        for per_ii=2:-1:1
            
             delta_edges=0.02;
            center_edges=[-0.5+(delta_edges/2):delta_edges:0.5-(delta_edges/2)];
            edges=[-0.5:delta_edges:0.5];
            Cxy_histo=zeros(200,length(center_edges));
            ii_Cxy_histo=0;
            
            grNo=1;
            
            bar_offset = bar_offset + 1;
            
            %Get these coherence values
         
            for ii=1:length(FileName)
                this_jj=[];
                for jj=1:all_files(ii).handles_out.dcoh_ii
                    
                    if all_files(ii).handles_out.odor_Cxy_values(jj).bwii==bwii
                        if all_files(ii).handles_out.odor_Cxy_values(jj).evNo==evNo
                            if all_files(ii).handles_out.odor_Cxy_values(jj).per_ii==per_ii
                                if all_files(ii).handles_out.odor_Cxy_values(jj).grNo==grNo
                                    this_jj=jj;
                                end
                            end
                            
                        end
                    end
                end
                if ~isempty(this_jj)
                    
                    for ii_mouse=1:length(all_files(ii).handles_out.odor_Cxy_values(this_jj).Cxy_per_mouse)
                        these_xCxy=all_files(ii).handles_out.odor_Cxy_values(this_jj).Cxy_per_mouse(ii_mouse).x_Cxy;
                        [N,edges]=histcounts(these_xCxy,edges);
                        if sum(isnan(N/sum(N)))==0
                            ii_Cxy_histo=ii_Cxy_histo+1;
                            Cxy_histo(ii_Cxy_histo,:)=N/sum(N);
                        end
                        %                         these_f_Cxy=all_files(ii).handles_out.odor_Cxy_values(this_jj).Cxy_per_mouse(ii_mouse).f_Cxy;
                        %                         last_f=0;
                        %                         ii_these_coh=ii_these_coh+1;
                        %                         for ii_c=1:length(coh_values)
                        %                             these_iis=find((these_xCxy>=coh_values(ii_c)-0.025)&(these_xCxy<coh_values(ii_c)+0.025));
                        %                             if isempty(these_iis)
                        %                                 these_coh(ii_these_coh,ii_c)=last_f;
                        %                             else
                        %                                 these_coh(ii_these_coh,ii_c)=mean(these_f_Cxy(these_iis));
                        %                             end
                        %                             last_f=these_coh(ii_these_coh,ii_c);
                        %                         end
                    end
                end
            end
            
            Cxy_histo=Cxy_histo(1:ii_Cxy_histo,:);
            
            mean_Cxy=[];
            mean_Cxy=mean(Cxy_histo,1);
            
            CI=[];
            CI = bootci(1000, {@mean, Cxy_histo})';
            maxlP=max([maxlP max(CI(:))]);
            minlP=min([minlP min(CI(:))]);
            CI(:,1)= mean_Cxy'-CI(:,1);
            CI(:,2)=CI(:,2)- mean_Cxy';
            
            
            
            
            if evNo==2
                if per_ii==1
                    %S- Proficient
                    [hlCR, hpCR] = boundedline(center_edges',mean_Cxy', CI, 'cmap',[158/255 31/255 99/255]);
                    %                     bar(bar_offset,mean(these_coh),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
                else
                    %S- Naive
                    [hlCR, hpCR] = boundedline(center_edges',mean_Cxy', CI, 'cmap',[238/255 111/255 179/255]);
                    %                     bar(bar_offset,mean(these_coh),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
                end
            else
                if per_ii==1
                    %S+ Proficient
                    [hlCR, hpCR] = boundedline(center_edges',mean_Cxy', CI, 'cmap',[0 114/255 178/255]);
                    %                     bar(bar_offset,mean(these_coh),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
                else
                    %S+ naive
                    [hlCR, hpCR] = boundedline(center_edges',mean_Cxy', CI, 'cmap',[80/255 194/255 255/255]);
                    %                     bar(bar_offset,mean(these_coh),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
                end
            end
            
            %             CI = bootci(1000, {@mean, these_coh},'type','cper');
            %             plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
            %             plot(bar_offset*ones(1,length(these_coh)),these_coh,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
            %
            
            %
            %             glm_coh.data(glm_ii+1:glm_ii+length(these_coh))=these_coh;
            %             %                 glm_coh.group(glm_ii+1:glm_ii+length(these_coh))=grNo*ones(1,length(these_coh));
            %             glm_coh.perCorr(glm_ii+1:glm_ii+length(these_coh))=per_ii*ones(1,length(these_coh));
            %             glm_coh.event(glm_ii+1:glm_ii+length(these_coh))=evNo*ones(1,length(these_coh));
            %             glm_ii=glm_ii+length(these_coh);
            %
            %             id_ii=id_ii+1;
            %             input_data(id_ii).data=these_coh;
            %             input_data(id_ii).description=[ evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
            
        end
        %         drgViolinPlot([bar_offset-1 bar_offset],these_coh_per_Ev,edges,rand_offset,'k','k',4,1);
        
        plot([0 0],[0 0.2],'-k')
        
        title(['Odor coherence for ' evTypeLabels{evNo} ' ' bandwidth_names{bwii}])
        
        
        
        xlim([-0.5 0.5])
        ylim([0 0.2])
        
        ylabel('Probability')
        xlabel('Coherence')
        
        
    end
    
    
    
    %     %Perform the glm
    %     fprintf(1, ['glm for odor coherence per mouse per odor pair for '  bandwidth_names{bwii} '\n'])
    %     fprintf(fileID, ['glm for odor coherence per mouse per odor pair for '  bandwidth_names{bwii} '\n']);
    %
    %     fprintf(1, ['\n\nglm for odor coherence for' bandwidth_names{bwii} '\n'])
    %     fprintf(fileID, ['\n\nglm for odor coherence for' bandwidth_names{bwii} '\n']);
    %
    %     tbl = table(glm_coh.data',glm_coh.perCorr',glm_coh.event',...
    %         'VariableNames',{'MI','perCorr','event'});
    %     mdl = fitglm(tbl,'MI~perCorr+event+perCorr*event'...
    %         ,'CategoricalVars',[2,3]);
    %
    %     txt = evalc('mdl');
    %     txt=regexp(txt,'<strong>','split');
    %     txt=cell2mat(txt);
    %     txt=regexp(txt,'</strong>','split');
    %     txt=cell2mat(txt);
    %
    %     fprintf(fileID,'%s\n', txt);
    %
    %
    %     %Do the ranksum/t-test
    %     fprintf(1, ['\n\nRanksum or t-test p values for odor coherence per mouse per odor pair for ' bandwidth_names{bwii} ' \n'])
    %     fprintf(fileID, ['\n\nRanksum or t-test p values for odor coherence per mouse per odor pair for ' bandwidth_names{bwii} ' \n']);
    %     [output_data] = drgMutiRanksumorTtest(input_data);
    
    
end

%Now plot the odor coherence per bandwidth

%Now plot the cumulative histos
edges=[-0.5:0.05:0.5];
rand_offset=0.8;


for bwii=[1 2 4]    %for amplitude bandwidths (beta, low gamma, high gamma)
    
    %Plot the average
    figNo = figNo +1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    ax=gca;ax.LineWidth=3;
    
    set(hFig, 'units','normalized','position',[.1 .5 .3 .4])
    hold on
    
    glm_coh=[];
    glm_ii=0;
    
    id_ii=0;
    input_data=[];
    
    
    
    ax=gca;ax.LineWidth=3;
    
    delta_edges=0.02;
    center_edges=[-0.5+(delta_edges/2):delta_edges:0.5-(delta_edges/2)];
    edges=[-0.5:delta_edges:0.5];
    Cxy_histo=zeros(200,length(center_edges));
    all_mean_Cxy=[];
    sh_Cxy_histo=zeros(200,length(center_edges));
    all_mean_shCxy=[];
    ii_Cxy_histo=0;
    ii_shCxy_histo=0;
    
    for evNo=1:2
        
        
        for per_ii=2:-1:1
            
            
            
            grNo=1;
            
            
            %Get these coherence values
            
            for ii=1:length(FileName)
                this_jj=[];
                for jj=1:all_files(ii).handles_out.dcoh_ii
                    
                    if all_files(ii).handles_out.odor_Cxy_values(jj).bwii==bwii
                        if all_files(ii).handles_out.odor_Cxy_values(jj).evNo==evNo
                            if all_files(ii).handles_out.odor_Cxy_values(jj).per_ii==per_ii
                                if all_files(ii).handles_out.odor_Cxy_values(jj).grNo==grNo
                                    this_jj=jj;
                                end
                            end
                            
                        end
                    end
                end
                if ~isempty(this_jj)
                    
                    for ii_mouse=1:length(all_files(ii).handles_out.odor_Cxy_values(this_jj).Cxy_per_mouse)
                        these_xCxy=all_files(ii).handles_out.odor_Cxy_values(this_jj).Cxy_per_mouse(ii_mouse).x_Cxy;
                        [N,edges]=histcounts(these_xCxy,edges);
                        if sum(isnan(N/sum(N)))==0
                            ii_Cxy_histo=ii_Cxy_histo+1;
                            Cxy_histo(ii_Cxy_histo,:)=N/sum(N);
                             all_mean_Cxy(ii_Cxy_histo)=mean(these_xCxy);
                        end
                        
                        these_sh_xCxy=all_files(ii).handles_out.odor_Cxy_values(this_jj).Cxy_per_mouse(ii_mouse).x_shCxy;
                        [N,edges]=histcounts(these_sh_xCxy,edges);
                        if sum(isnan(N/sum(N)))==0
                            ii_shCxy_histo=ii_shCxy_histo+1;
                            shCxy_histo(ii_shCxy_histo,:)=N/sum(N);
                            all_mean_shCxy(ii_shCxy_histo)=mean(these_sh_xCxy);
                        end
                        
                    end
                end
            end
            
            
            %             CI = bootci(1000, {@mean, these_coh},'type','cper');
            %             plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
            %             plot(bar_offset*ones(1,length(these_coh)),these_coh,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
            %
            
            %
            %             glm_coh.data(glm_ii+1:glm_ii+length(these_coh))=these_coh;
            %             %                 glm_coh.group(glm_ii+1:glm_ii+length(these_coh))=grNo*ones(1,length(these_coh));
            %             glm_coh.perCorr(glm_ii+1:glm_ii+length(these_coh))=per_ii*ones(1,length(these_coh));
            %             glm_coh.event(glm_ii+1:glm_ii+length(these_coh))=evNo*ones(1,length(these_coh));
            %             glm_ii=glm_ii+length(these_coh);
            %
            %             id_ii=id_ii+1;
            %             input_data(id_ii).data=these_coh;
            %             input_data(id_ii).description=[ evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
            
        end
        %         drgViolinPlot([bar_offset-1 bar_offset],these_coh_per_Ev,edges,rand_offset,'k','k',4,1);
        
        
        
    end
    
    
    
    %     %Perform the glm
    %     fprintf(1, ['glm for odor coherence per mouse per odor pair for '  bandwidth_names{bwii} '\n'])
    %     fprintf(fileID, ['glm for odor coherence per mouse per odor pair for '  bandwidth_names{bwii} '\n']);
    %
    %     fprintf(1, ['\n\nglm for odor coherence for' bandwidth_names{bwii} '\n'])
    %     fprintf(fileID, ['\n\nglm for odor coherence for' bandwidth_names{bwii} '\n']);
    %
    %     tbl = table(glm_coh.data',glm_coh.perCorr',glm_coh.event',...
    %         'VariableNames',{'MI','perCorr','event'});
    %     mdl = fitglm(tbl,'MI~perCorr+event+perCorr*event'...
    %         ,'CategoricalVars',[2,3]);
    %
    %     txt = evalc('mdl');
    %     txt=regexp(txt,'<strong>','split');
    %     txt=cell2mat(txt);
    %     txt=regexp(txt,'</strong>','split');
    %     txt=cell2mat(txt);
    %
    %     fprintf(fileID,'%s\n', txt);
    %
    %
    %     %Do the ranksum/t-test
    %     fprintf(1, ['\n\nRanksum or t-test p values for odor coherence per mouse per odor pair for ' bandwidth_names{bwii} ' \n'])
    %     fprintf(fileID, ['\n\nRanksum or t-test p values for odor coherence per mouse per odor pair for ' bandwidth_names{bwii} ' \n']);
    %     [output_data] = drgMutiRanksumorTtest(input_data);
    shCxy_histo=shCxy_histo(1:ii_shCxy_histo,:);
    
    mean_Cxy=[];
    mean_Cxy=mean(shCxy_histo,1);
    max_mean_Cxy=max(mean_Cxy);
    mean_Cxy=mean_Cxy/max_mean_Cxy;
    shCxy_histo=shCxy_histo/max_mean_Cxy;
    
    CI=[];
    CI = bootci(1000, {@mean, shCxy_histo})';
    maxlP=max([maxlP max(CI(:))]);
    minlP=min([minlP min(CI(:))]);
    CI(:,1)= mean_Cxy'-CI(:,1);
    CI(:,2)=CI(:,2)- mean_Cxy';
    
    [hlCR, hpCR] = boundedline(center_edges',mean_Cxy', CI, 'k');
    
    sh_mean_Cxy=mean_Cxy;
    
    plot([mean(all_mean_shCxy) mean(all_mean_shCxy)],[0 1.1],'-k','LineWidth',2)
    
    Cxy_histo=Cxy_histo(1:ii_Cxy_histo,:);
    
    mean_Cxy=[];
    mean_Cxy=mean(Cxy_histo,1);
    max_mean_Cxy=max(mean_Cxy);
    mean_Cxy=mean_Cxy/max_mean_Cxy;
    Cxy_histo=Cxy_histo/max_mean_Cxy;
    
    CI=[];
    CI = bootci(1000, {@mean, Cxy_histo})';
    maxlP=max([maxlP max(CI(:))]);
    minlP=min([minlP min(CI(:))]);
    CI(:,1)= mean_Cxy'-CI(:,1);
    CI(:,2)=CI(:,2)- mean_Cxy';
    
    [hlCR, hpCR] = boundedline(center_edges',mean_Cxy', CI, 'm');
    
    plot([mean(all_mean_Cxy) mean(all_mean_Cxy)],[0 1.1],'-m','LineWidth',2)
    
    plot(center_edges,sh_mean_Cxy,'-k','LineWidth',3)
    plot(center_edges,mean_Cxy,'-m','LineWidth',3)
    
    title(['Odor coherence for ' bandwidth_names{bwii}])
    
    xlim([-0.5 0.5])
    ylim([0 1.2])
    
    ylabel('A.U.')
    xlabel('Coherence')
    
end



%Now plot the percent significant delta Cxy
edges=[0:5:100];
rand_offset=0.8;


for bwii=[1 2 4]    %for amplitude bandwidths (beta, low gamma, high gamma)
    
    glm_coh=[];
    glm_ii=0;
    
    id_ii=0;
    input_data=[];
    
     glm_coh_mm=[];
    glm_ii_mm=0;
    
    id_ii_mm=0;
    input_data_mm=[];
    
    %Plot the average
    figNo = figNo +1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    ax=gca;ax.LineWidth=3;
    
    set(hFig, 'units','normalized','position',[.1 .5 .3 .4])
    hold on
    
    ax=gca;ax.LineWidth=3;
    
    bar_lab_loc=[];
    no_ev_labels=0;
    ii_gr_included=0;
    bar_offset = 0;
    
    for evNo=1:2
        
%         these_coh_per_Ev=[];
%         ii_coh=0;
        per_sig_thisEv=[];
        mouse_no_thisEv=[];
        per_ii_thisEv=[];
        
        for per_ii=2:-1:1
            
            grNo=1;
            
            bar_offset = bar_offset + 1;
            
            %Get these coherence values
            these_per_sig=[];
            ii_these_per_sig=0;
            for ii=1:length(FileName)
                these_mouse_no_per_op=T_mouse_no.mouse_no_per_op(T_mouse_no.odor_pair_no==ii);
                these_mouse_no=T_mouse_no.mouse_no(T_mouse_no.odor_pair_no==ii);
                this_jj=[];
                for jj=1:all_files(ii).handles_out.dcoh_ii
                    
                    if all_files(ii).handles_out.per_sig_values(jj).bwii==bwii
                        if all_files(ii).handles_out.per_sig_values(jj).evNo==evNo
                            if all_files(ii).handles_out.per_sig_values(jj).per_ii==per_ii
                                
                                
                                if all_files(ii).handles_out.per_sig_values(jj).grNo==grNo
                                    this_jj=jj;
                                end
                                
                            end
                            
                        end
                    end
                end
                if ~isempty(this_jj)
                    if mouse_op==1
                        these_per_sig(ii_these_per_sig+1:ii_these_per_sig+length(all_files(ii).handles_out.per_sig_values(this_jj).percent_sig_per_mouse))=all_files(ii).handles_out.per_sig_values(this_jj).percent_sig_per_mouse;
                        ii_these_per_sig=ii_these_per_sig+length(all_files(ii).handles_out.per_sig_values(this_jj).percent_sig_per_mouse);
                    else
                        ii_these_per_sig=ii_these_per_sig+1;
                        these_per_sig(ii_these_per_sig)=all_files(ii).handles_out.per_sig_values(this_jj).dcoh;
                    end
                    
                    %Assign to overall mouse numbers
                    these_percent_sig=all_files(ii).handles_out.per_sig_values(this_jj).percent_sig_per_mouse;
                    these_mouse_nos=all_files(ii).handles_out.per_sig_values(this_jj).mouseNo;
                    these_mouse_no_overall=[];
                    for ii_mouse=1:length(these_mouse_nos)
                        these_mouse_no_overall(ii_mouse)=these_mouse_no(find(these_mouse_nos(ii_mouse)==these_mouse_no_per_op));
                    end
                    per_sig_thisEv=[per_sig_thisEv these_percent_sig];
                    mouse_no_thisEv=[mouse_no_thisEv these_mouse_no_overall];
                    per_ii_thisEv=[per_ii_thisEv per_ii*ones(1,length(these_percent_sig))];

                end
            end
            
              
            if evNo==2
                if per_ii==1
                    %S- Proficient
                    bar(bar_offset,mean(these_per_sig),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
%                     drgViolinPlot([bar_offset-1 bar_offset],[last_these_per_sig;these_per_sig],edges,rand_offset,'k','k',4,1);
                else
                    %S- Naive
                    bar(bar_offset,mean(these_per_sig),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
                    last_these_per_sig=these_per_sig;
                end
            else
                if per_ii==1
                    %S+ Proficient
                    bar(bar_offset,mean(these_per_sig),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
%                     drgViolinPlot([bar_offset-1 bar_offset],[last_these_per_sig;these_per_sig],edges,rand_offset,'k','k',4,1);
                else
                    %S+ naive
                    bar(bar_offset,mean(these_per_sig),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
                    last_these_per_sig=these_per_sig;
                end
            end
            
%             ii_these_per_sig=ii_these_per_sig+1;
%             these_per_sig_per_Ev(ii_these_per_sig,:)=these_per_sig;
            
            %Violin plot
            
%             [mean_out, CIout]=drgViolinPoint(these_per_sig,edges,bar_offset,rand_offset,'k','k',3);
            %             CI = bootci(1000, {@mean, these_per_sig},'type','cper');
            %             plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
            %             plot(bar_offset*ones(1,length(these_per_sig)),these_per_sig,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
            %
            
            
            glm_coh.data(glm_ii+1:glm_ii+length(these_per_sig))=these_per_sig;
            glm_coh.mouse_no(glm_ii+1:glm_ii+length(these_per_sig))=mouse_no_thisEv(per_ii_thisEv==per_ii);
            glm_coh.perCorr(glm_ii+1:glm_ii+length(these_per_sig))=per_ii*ones(1,length(these_per_sig));
            glm_coh.event(glm_ii+1:glm_ii+length(these_per_sig))=evNo*ones(1,length(these_per_sig));
            glm_ii=glm_ii+length(these_per_sig);
            
            id_ii=id_ii+1;
            input_data(id_ii).data=these_per_sig;
            input_data(id_ii).description=[ evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
             
        end
%         drgViolinPlot([bar_offset-1 bar_offset],these_per_sig_per_Ev,edges,rand_offset,'k','k',4,1);
        %Calculate the per mouse means

        per_ii2_mean_per_sig_mm=[];
        per_ii1_mean_per_sig_mm=[];
        for mouseNo=1:max(mouse_no_thisEv)
            per_ii1_mean_per_sig_mm=[per_ii1_mean_per_sig_mm mean(per_sig_thisEv((mouse_no_thisEv==mouseNo)&(per_ii_thisEv==1)))];
            per_ii2_mean_per_sig_mm=[per_ii2_mean_per_sig_mm mean(per_sig_thisEv((mouse_no_thisEv==mouseNo)&(per_ii_thisEv==2)))];
        end
        
        %per_ii=1
        glm_coh_mm.data(glm_ii_mm+1:glm_ii_mm+length(per_ii1_mean_per_sig_mm))=per_ii1_mean_per_sig_mm;
        %                 glm_coh_mm.group(glm_ii_mm+1:glm_ii_mm+length(these_per_sig))=grNo*ones(1,length(per_ii1_mean_per_sig_mm));
        glm_coh_mm.perCorr(glm_ii_mm+1:glm_ii_mm+length(per_ii1_mean_per_sig_mm))=ones(1,length(per_ii1_mean_per_sig_mm));
        glm_coh_mm.event(glm_ii_mm+1:glm_ii_mm+length(per_ii1_mean_per_sig_mm))=evNo*ones(1,length(per_ii1_mean_per_sig_mm));
        glm_ii_mm=glm_ii_mm+length(per_ii1_mean_per_sig_mm);
        
        id_ii_mm=id_ii_mm+1;
        input_data_mm(id_ii_mm).data=per_ii1_mean_per_sig_mm;
        input_data_mm(id_ii_mm).description=[ evTypeLabels{evNo} ' ' prof_naive_leg{1}];
        
         %per_ii=1
        glm_coh_mm.data(glm_ii_mm+1:glm_ii_mm+length(per_ii2_mean_per_sig_mm))=per_ii2_mean_per_sig_mm;
        %                 glm_coh_mm.group(glm_ii_mm+1:glm_ii_mm+length(these_per_sig))=grNo*ones(1,length(per_ii2_mean_per_sig_mm));
        glm_coh_mm.perCorr(glm_ii_mm+1:glm_ii_mm+length(per_ii2_mean_per_sig_mm))=2*ones(1,length(per_ii2_mean_per_sig_mm));
        glm_coh_mm.event(glm_ii_mm+1:glm_ii_mm+length(per_ii2_mean_per_sig_mm))=evNo*ones(1,length(per_ii2_mean_per_sig_mm));
        glm_ii_mm=glm_ii_mm+length(per_ii2_mean_per_sig_mm);
        
        id_ii_mm=id_ii_mm+1;
        input_data_mm(id_ii_mm).data=per_ii2_mean_per_sig_mm;
        input_data_mm(id_ii_mm).description=[ evTypeLabels{evNo} ' ' prof_naive_leg{2}];
        
        for mouseNo=1:length(per_ii2_mean_per_sig_mm)
            plot([bar_offset-1 bar_offset],[per_ii2_mean_per_sig_mm(mouseNo) per_ii1_mean_per_sig_mm(mouseNo)],'-ok','MarkerSize',6,'LineWidth',2,'MarkerFaceColor',[0.6 0.6 0.6],'MarkerEdgeColor',[0.6 0.6 0.6],'Color',[0.6 0.6 0.6])
        end
        
        %Violin plot
        [mean_out, CIout]=drgViolinPoint(these_per_sig...
            ,edges,bar_offset,rand_offset,'k','k',3);
        
        %Violin plot
        [mean_out, CIout]=drgViolinPoint(last_these_per_sig...
            ,edges,bar_offset-1,rand_offset,'k','k',3);
        
        bar_offset = bar_offset + 1;
        
    end
    
    title(['Percent significant delta coherence for ' bandwidth_names{bwii}])
    ylim([0 100])
    
    %     %Annotations identifying groups
    %     x_interval=0.8/ii_gr_included;
    %     for ii=1:ii_gr_included
    %         annotation('textbox',[0.7*x_interval+x_interval*(ii-1) 0.7 0.3 0.1],'String',handles_drgb.drgbchoices.group_no_names{ groups_included(ii)},'FitBoxToText','on');
    %     end
    %
    %     %Proficient/Naive annotations
    %     annotation('textbox',[0.15 0.8 0.3 0.1],'String','Proficient','FitBoxToText','on','Color','r','LineStyle','none');
    %     annotation('textbox',[0.15 0.75 0.3 0.1],'String','Naive','FitBoxToText','on','Color','b','LineStyle','none');
    
    plot([0 6],[0 0],'-k','LineWidth',2)
    xlim([0 6])
    
    xticks([1 2 4 5])
    xticklabels({'nS+', 'pS+','nS-', 'pS-'})
    
    ylabel('Percent significant')
    
    
    %Perform the glm per mouse per odor pair
    fprintf(1, ['glm for percent significant delta imaginary coherence per mouse per odor pair for '  bandwidth_names{bwii} '\n'])
    fprintf(fileID, ['glm for percent significant delta imaginary coherence per mouse per odor pair for '  bandwidth_names{bwii} '\n']);
    
    tbl = table(glm_coh.data',glm_coh.perCorr',glm_coh.event',...
        'VariableNames',{'per_sig','naive_vs_proficient','sp_vs_sm'});
    mdl = fitglm(tbl,'per_sig~naive_vs_proficient+sp_vs_sm+naive_vs_proficient*sp_vs_sm'...
        ,'CategoricalVars',[2,3])
    
    txt = evalc('mdl');
    txt=regexp(txt,'<strong>','split');
    txt=cell2mat(txt);
    txt=regexp(txt,'</strong>','split');
    txt=cell2mat(txt);
    
    fprintf(fileID,'%s\n', txt);
    
    
    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for percent significant delta imaginary coherence per mouse per odor pair for ' bandwidth_names{bwii} ' \n'])
    fprintf(fileID, ['\n\nRanksum or t-test p values for percent significant delta imaginary coherence per mouse per odor pair for ' bandwidth_names{bwii} ' \n']);
    [output_data] = drgMutiRanksumorTtest(input_data, fileID);
    
    %Nested ANOVAN
    %https://www.mathworks.com/matlabcentral/answers/491365-within-between-subjects-in-anovan
    nesting=[0 0 0; ... % This line indicates that group factor is not nested in any other factor.
         0 0 0; ... % This line indicates that perCorr is not nested in any other factor.
         1 1 0];    % This line indicates that mouse_no (the last factor) is nested under group, perCorr and event
                    % (the 1 in position 1 on the line indicates nesting under the first factor).
    figNo=figNo+1;
                     
    [p anovanTbl stats]=anovan(glm_coh.data,{glm_coh.perCorr glm_coh.event glm_coh.mouse_no},...
    'model','interaction',...
    'nested',nesting,...
    'varnames',{'naive_vs_proficient', 'sp_vs_sm','mouse_no'});
  
    fprintf(fileID, ['\n\nNested ANOVAN for percent significant delta imaginary coherence per mouse per odor pair for '  bandwidth_names{bwii} '\n'])
    drgWriteANOVANtbl(anovanTbl,fileID);
    fprintf(fileID, '\n\n')
    
    
    %Perform the glm per mouse 
    fprintf(1, ['glm for percent significant delta imaginary coherence per mouse for '  bandwidth_names{bwii} '\n'])
%     fprintf(fileID, ['glm for percent significant delta imaginary coherence per mouse for '  bandwidth_names{bwii} '\n']);
    
    tbl = table(glm_coh_mm.data',glm_coh_mm.perCorr',glm_coh_mm.event',...
        'VariableNames',{'per_sig','naive_vs_proficient','sp_vs_sm'});
    mdl = fitglm(tbl,'per_sig~naive_vs_proficient+sp_vs_sm+naive_vs_proficient*sp_vs_sm'...
        ,'CategoricalVars',[2,3])
    
%     txt = evalc('mdl');
%     txt=regexp(txt,'<strong>','split');
%     txt=cell2mat(txt);
%     txt=regexp(txt,'</strong>','split');
%     txt=cell2mat(txt);
%     
%     fprintf(fileID,'%s\n', txt);
    
    
    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for percent significant delta imaginary coherence per mouse for ' bandwidth_names{bwii} ' \n'])
%     fprintf(fileID, ['\n\nRanksum or t-test p values for percent significant delta imaginary coherence per mouse for ' bandwidth_names{bwii} ' \n']);
    [output_data] = drgMutiRanksumorTtest(input_data_mm);
    
end

%Now plot the deltaCxy histogram for significant deltaCxy

rand_offset=0.8;


for bwii=[1 2 4]    %for amplitude bandwidths (beta, low gamma, high gamma)
    
    glm_coh=[];
    glm_ii=0;
    
    id_ii=0;
    input_data=[];
    
    %Plot the average
    figNo = figNo +1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    ax=gca;ax.LineWidth=3;
    
    set(hFig, 'units','normalized','position',[.1 .5 .3 .4])
    hold on
    
    ax=gca;ax.LineWidth=3;
    
    bar_lab_loc=[];
    no_ev_labels=0;
    ii_gr_included=0;
    bar_offset = 0;
    
    for evNo=1:2
        
        these_coh_per_Ev=[];
        ii_coh=0;
        
        for per_ii=2:-1:1
            
            delta_edges=0.02;
            center_edges=[-0.5+(delta_edges/2):delta_edges:0.5-(delta_edges/2)];
            edges=[-0.5:delta_edges:0.5];
            deltaCxy_histo=zeros(200,length(center_edges));
            ii_deltaCxy_histo=0;
            grNo=1;
            
            bar_offset = bar_offset + 1;
            
            %Get these deltaCxy values
            for ii=1:length(FileName)
                this_jj=[];
                for jj=1:all_files(ii).handles_out.dcoh_ii
                    
                    if all_files(ii).handles_out.per_sig_values(jj).bwii==bwii
                        if all_files(ii).handles_out.per_sig_values(jj).evNo==evNo
                            if all_files(ii).handles_out.per_sig_values(jj).per_ii==per_ii
                                if all_files(ii).handles_out.per_sig_values(jj).grNo==grNo
                                    this_jj=jj;
                                end
                            end
                        end
                    end
                end
                if ~isempty(this_jj)
                    if mouse_op==1
                        for ww=1:length(all_files(ii).handles_out.per_sig_values(this_jj).sig_deltaCxy_per_mouse)
                            [N,edges]=histcounts(all_files(ii).handles_out.per_sig_values(this_jj).sig_deltaCxy_per_mouse(ww).deltaCxy,edges);
                            if sum(isnan(N/sum(N)))==0
                                ii_deltaCxy_histo=ii_deltaCxy_histo+1;
                                deltaCxy_histo(ii_deltaCxy_histo,:)=N/sum(N);
                            end
                        end
                    end
                    
                end
            end
            
            deltaCxy_histo=deltaCxy_histo(1:ii_deltaCxy_histo,:);
            
            mean_deltaCxy=[];
            mean_deltaCxy=mean(deltaCxy_histo,1);
            
            CI=[];
            CI = bootci(1000, {@mean, deltaCxy_histo})';
            maxlP=max([maxlP max(CI(:))]);
            minlP=min([minlP min(CI(:))]);
            CI(:,1)= mean_deltaCxy'-CI(:,1);
            CI(:,2)=CI(:,2)- mean_deltaCxy';
            
            
            if evNo==2
                if per_ii==1
                    %S- Proficient
                    [hlCR, hpCR] = boundedline(center_edges',mean_deltaCxy', CI, 'cmap',[158/255 31/255 99/255]);
                    %                     bar(bar_offset,mean(these_coh),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
                else
                    %S- Naive
                    [hlCR, hpCR] = boundedline(center_edges',mean_deltaCxy', CI, 'cmap',[238/255 111/255 179/255]);
                    %                     bar(bar_offset,mean(these_coh),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
                end
            else
                if per_ii==1
                    %S+ Proficient
                    [hlCR, hpCR] = boundedline(center_edges',mean_deltaCxy', CI, 'cmap',[0 114/255 178/255]);
                    %                     bar(bar_offset,mean(these_coh),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
                else
                    %S+ naive
                    [hlCR, hpCR] = boundedline(center_edges',mean_deltaCxy', CI, 'cmap',[80/255 194/255 255/255]);
                    %                     bar(bar_offset,mean(these_coh),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
                end
            end
            
            
            
            %             glm_coh.data(glm_ii+1:glm_ii+length(these_per_sig))=these_per_sig;
            %             %                 glm_coh.group(glm_ii+1:glm_ii+length(these_per_sig))=grNo*ones(1,length(these_per_sig));
            %             glm_coh.perCorr(glm_ii+1:glm_ii+length(these_per_sig))=per_ii*ones(1,length(these_per_sig));
            %             glm_coh.event(glm_ii+1:glm_ii+length(these_per_sig))=evNo*ones(1,length(these_per_sig));
            %             glm_ii=glm_ii+length(these_per_sig);
            %
            %             id_ii=id_ii+1;
            %             input_data(id_ii).data=these_per_sig;
            %             input_data(id_ii).description=[ evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
            
        end
        
    end
    
    title(['Histogram for delta odor coherence for significant changes for ' bandwidth_names{bwii}])
    
    xlabel('Delta coherence')
    ylabel('Fraction')
    
    
    %     %Perform the glm
    %     fprintf(1, ['glm for odor coherence per mouse per odor pair for '  bandwidth_names{bwii} '\n'])
    %     fprintf(fileID, ['glm for odor coherence per mouse per odor pair for '  bandwidth_names{bwii} '\n']);
    %
    %     fprintf(1, ['\n\nglm for odor coherence for' bandwidth_names{bwii} '\n'])
    %     fprintf(fileID, ['\n\nglm for odor coherence for' bandwidth_names{bwii} '\n']);
    %
    %     tbl = table(glm_coh.data',glm_coh.perCorr',glm_coh.event',...
    %         'VariableNames',{'MI','perCorr','event'});
    %     mdl = fitglm(tbl,'MI~perCorr+event+perCorr*event'...
    %         ,'CategoricalVars',[2,3])
    %
    %     txt = evalc('mdl');
%     txt=regexp(txt,'<strong>','split');
%     txt=cell2mat(txt);
%     txt=regexp(txt,'</strong>','split');
%     txt=cell2mat(txt);
%     
%     fprintf(fileID,'%s\n', txt);
%     
%     
%     %Do the ranksum/t-test
%     fprintf(1, ['\n\nRanksum or t-test p values for odor coherence per mouse per odor pair for ' bandwidth_names{bwii} ' \n'])
%     fprintf(fileID, ['\n\nRanksum or t-test p values for odor coherence per mouse per odor pair for ' bandwidth_names{bwii} ' \n']);
%     [output_data] = drgMutiRanksumorTtest(input_data);
    
    
end

%Now plot the histograms for odor Cxy
rand_offset=0.8;


for bwii=[1 2 4]    %for amplitude bandwidths (beta, low gamma, high gamma)
    
    glm_coh=[];
    glm_ii=0;
    
    id_ii=0;
    input_data=[];
    
    %Plot the average
    figNo = figNo +1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    ax=gca;ax.LineWidth=3;
    
    set(hFig, 'units','normalized','position',[.1 .5 .3 .4])
    hold on
    
    ax=gca;ax.LineWidth=3;
    
    bar_lab_loc=[];
    no_ev_labels=0;
    ii_gr_included=0;
    bar_offset = 0;
    
    for evNo=1:2
        
        these_coh_per_Ev=[];
        ii_coh=0;
        
        for per_ii=2:-1:1
            
            delta_edges=0.02;
            center_edges=[-0.5+(delta_edges/2):delta_edges:0.5-(delta_edges/2)];
            edges=[-0.5:delta_edges:0.5];
            Cxy_histo=zeros(200,length(center_edges));
            ii_Cxy_histo=0;
            grNo=1;
            
            bar_offset = bar_offset + 1;
            
            %Get these Cxy values
            for ii=1:length(FileName)
                this_jj=[];
                for jj=1:all_files(ii).handles_out.dcoh_ii
                    
                    if all_files(ii).handles_out.per_sig_values(jj).bwii==bwii
                        if all_files(ii).handles_out.per_sig_values(jj).evNo==evNo
                            if all_files(ii).handles_out.per_sig_values(jj).per_ii==per_ii
                                if all_files(ii).handles_out.per_sig_values(jj).grNo==grNo
                                    this_jj=jj;
                                end
                            end
                        end
                    end
                end
                if ~isempty(this_jj)
                    if mouse_op==1
                        for ww=1:length(all_files(ii).handles_out.per_sig_values(this_jj).sig_Cxy_per_mouse)
                            [N,edges]=histcounts(all_files(ii).handles_out.per_sig_values(this_jj).sig_Cxy_per_mouse(ww).Cxy,edges);
                            if sum(isnan(N/sum(N)))==0
                                ii_Cxy_histo=ii_Cxy_histo+1;
                                Cxy_histo(ii_Cxy_histo,:)=N/sum(N);
                            end
                        end
                    end
                    
                end
            end
            
            Cxy_histo=Cxy_histo(1:ii_Cxy_histo,:);
            
            mean_Cxy=[];
            mean_Cxy=mean(Cxy_histo,1);
            
            CI=[];
            CI = bootci(1000, {@mean, Cxy_histo})';
            maxlP=max([maxlP max(CI(:))]);
            minlP=min([minlP min(CI(:))]);
            CI(:,1)= mean_Cxy'-CI(:,1);
            CI(:,2)=CI(:,2)- mean_Cxy';
            
            
            if evNo==2
                if per_ii==1
                    %S- Proficient
                    [hlCR, hpCR] = boundedline(center_edges',mean_Cxy', CI, 'cmap',[158/255 31/255 99/255]);
                    %                     bar(bar_offset,mean(these_coh),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
                else
                    %S- Naive
                    [hlCR, hpCR] = boundedline(center_edges',mean_Cxy', CI, 'cmap',[238/255 111/255 179/255]);
                    %                     bar(bar_offset,mean(these_coh),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
                end
            else
                if per_ii==1
                    %S+ Proficient
                    [hlCR, hpCR] = boundedline(center_edges',mean_Cxy', CI, 'cmap',[0 114/255 178/255]);
                    %                     bar(bar_offset,mean(these_coh),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
                else
                    %S+ naive
                    [hlCR, hpCR] = boundedline(center_edges',mean_Cxy', CI, 'cmap',[80/255 194/255 255/255]);
                    %                     bar(bar_offset,mean(these_coh),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
                end
            end
            
            
            
            %             glm_coh.data(glm_ii+1:glm_ii+length(these_per_sig))=these_per_sig;
            %             %                 glm_coh.group(glm_ii+1:glm_ii+length(these_per_sig))=grNo*ones(1,length(these_per_sig));
            %             glm_coh.perCorr(glm_ii+1:glm_ii+length(these_per_sig))=per_ii*ones(1,length(these_per_sig));
            %             glm_coh.event(glm_ii+1:glm_ii+length(these_per_sig))=evNo*ones(1,length(these_per_sig));
            %             glm_ii=glm_ii+length(these_per_sig);
            %
            %             id_ii=id_ii+1;
            %             input_data(id_ii).data=these_per_sig;
            %             input_data(id_ii).description=[ evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
            
        end
        
    end
    
    title(['Histogram for odor coherence for significant changes for ' bandwidth_names{bwii}])
    
    xlabel('Coherence')
    ylabel('Fraction')
    
    
    %     %Perform the glm
    %     fprintf(1, ['glm for odor coherence per mouse per odor pair for '  bandwidth_names{bwii} '\n'])
    %     fprintf(fileID, ['glm for odor coherence per mouse per odor pair for '  bandwidth_names{bwii} '\n']);
    %
    %     fprintf(1, ['\n\nglm for odor coherence for' bandwidth_names{bwii} '\n'])
    %     fprintf(fileID, ['\n\nglm for odor coherence for' bandwidth_names{bwii} '\n']);
    %
    %     tbl = table(glm_coh.data',glm_coh.perCorr',glm_coh.event',...
    %         'VariableNames',{'MI','perCorr','event'});
    %     mdl = fitglm(tbl,'MI~perCorr+event+perCorr*event'...
    %         ,'CategoricalVars',[2,3])
    %
    %     txt = evalc('mdl');
%     txt=regexp(txt,'<strong>','split');
%     txt=cell2mat(txt);
%     txt=regexp(txt,'</strong>','split');
%     txt=cell2mat(txt);
%     
%     fprintf(fileID,'%s\n', txt);
%     
%     
%     %Do the ranksum/t-test
%     fprintf(1, ['\n\nRanksum or t-test p values for odor coherence per mouse per odor pair for ' bandwidth_names{bwii} ' \n'])
%     fprintf(fileID, ['\n\nRanksum or t-test p values for odor coherence per mouse per odor pair for ' bandwidth_names{bwii} ' \n']);
%     [output_data] = drgMutiRanksumorTtest(input_data);
    
    
end
fclose(fileID);

pffft=1;


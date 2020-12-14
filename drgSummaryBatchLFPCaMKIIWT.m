function drgSummaryBatchLFPCaMKIIWT
%Analyzes the linear discriminant analysis performed by drgLFPDiscriminantBatch
%Takes as a in input the 'drgbDiscPar' file listing 'Discriminant_*.mat' output files
%from drgLFPDiscriminantBatch
%
%Performs summary analises for LDA and  PAC analysis performed with case 19
%of drgAnalysisBatchLFP

warning('off')

close all
clear all

frequency=[4:(95-4)/100:95];

PACnames{1}='Beta';
PACnames{2}='Low gamma';
PACnames{3}='High gamma';

prof_naive_leg{1}='Proficient';
prof_naive_leg{2}='Naive';

group_legend{1}='WT';
group_legend{2}='Het';
group_legend{3}='KO';

evTypeLabels{1}='S+';
evTypeLabels{2}='S-';

%Location of files
PathName='F:\Datos summary CaMKII111720\LFPPAC drgAnalysisBatchLFPCaMKII case 20 output\';

%Hippocampus
hippFileName{1}='CaMKIIacetoLFPPACall121319_hippoLFP.mat';
hippFileName{2}='CaMKIIethylbenLFPPACal222019_hippoLFP.mat';
hippFileName{3}='CaMKIIethylaceLFPPACall12219_hippoLFP.mat';
hippFileName{4}='CaMKIIpropylaceLFPPACall1220_hippoLFP.mat';
hippFileName{5}='CaMKIIpointzero1ethylacepropylaceLFPPACall122019_hippoLFP.mat';
hippFileName{6}='CaMKIIpointzero1propylaceLFPPACall11620_hippoLFP.mat';
hippFileName{7}='CaMKIIpzz1ethylaceLFPPACall121119_hippoLFP.mat';
hippFileName{8}='CaMKIIpzz1propylaceLFPPACall121919_hippoLFP.mat';

%Prefrontal
preFileName{1}='CaMKIIacetoLFPPACall121319_prefronLFP.mat';
preFileName{2}='CaMKIIethylbenLFPPACal222019_prefronLFP.mat';
preFileName{3}='CaMKIIethylaceLFPPACall12219_prefronLFP.mat';
preFileName{4}='CaMKIIpropylaceLFPPACall1220_prefronLFP.mat';
preFileName{5}='CaMKIIpointzero1ethylacepropylaceLFPPACall122019_prefronLFP.mat';
preFileName{6}='CaMKIIpointzero1propylaceLFPPACall11620_prefronLFP.mat';
preFileName{7}='CaMKIIpzz1ethylaceLFPPACall121119_prefronLFP.mat';
preFileName{8}='CaMKIIpzz1propylaceLFPPACall121919_prefronLFP.mat';





%Load data hippocampus
all_hippo=[];

for ii=1:length(hippFileName)
    load([PathName hippFileName{ii}])
    all_hippo(ii).handles_out=handles_out;
end

%Load data prefrontal
all_pre=[];

for ii=1:length(preFileName)
    load([PathName preFileName{ii}])
    all_pre(ii).handles_out=handles_out;
end

figNo=0;





%Now plot for the hippocampus the dependence of delta LFP to frequency

glm_dB_hipp=[];
glm_ii_hipp=0;

glm_dB_pre=[];
glm_ii_pre=0;

glm_dB_both=[];
glm_ii_both=0;

maxdB=-1000000;
mindB=1000000;

for grNo=1:3
    %Plot the average
    figNo = figNo +1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    
    
    set(hFig, 'units','normalized','position',[.1 .5 .35 .3])
    hold on
    
    
    for evNo=1:2
        
        for per_ii=2:-1:1
            
            
            
            %Get these MI values
            these_dBWB=[];
            ii_dB=0;
            for ii=1:length(hippFileName)
                this_jj=[];
                for jj=1:length(all_hippo(ii).handles_out.delta_dB_WB)
                    if all_hippo(ii).handles_out.delta_dB_WB(jj).per_ii==per_ii
                        if all_hippo(ii).handles_out.delta_dB_WB(jj).evNo==evNo
                            if all_hippo(ii).handles_out.delta_dB_WB(jj).grNo==grNo
                                this_jj=jj;
                            end
                        end
                    end
                end
                if ~isempty(this_jj)
                    ii_dB=ii_dB+1;
                    these_dBWB(ii_dB,1:length(frequency))=all_hippo(ii).handles_out.delta_dB_WB(this_jj).mean_dbWB;
                end
            end
            
            this_mean_dbWB=zeros(1,length(frequency));
            this_mean_dbWB(1,:)=mean(these_dBWB,1);
            
            CI=[];
            CI = bootci(1000, {@mean, these_dBWB})';
            CI(:,1)= this_mean_dbWB'-CI(:,1);
            CI(:,2)=CI(:,2)- this_mean_dbWB';
            
            if evNo==1
                if per_ii==1
                    %S+ Proficient
                    [hlCR, hpCR] = boundedline(frequency',this_mean_dbWB', CI, 'cmap',[158/255 31/255 99/255]);
                else
                    %S+ Naive
                    [hlCR, hpCR] = boundedline(frequency',this_mean_dbWB', CI, 'cmap',[238/255 111/255 179/255]);
                end
            else
                if per_ii==1
                    %S- Proficient
                    [hlCR, hpCR] = boundedline(frequency',this_mean_dbWB', CI, 'cmap',[0 114/255 178/255]);
                else
                    %S- naive
                    [hlCR, hpCR] = boundedline(frequency',this_mean_dbWB', CI, 'cmap',[80/255 194/255 255/255]);
                end
            end
            
            no_dbWB=size(these_dBWB,1);
            for iif=1:length(frequency)
                glm_dB_hipp.data(glm_ii_hipp+1:glm_ii_hipp+no_dbWB)=these_dBWB(:,iif);
                glm_dB_hipp.freq(glm_ii_hipp+1:glm_ii_hipp+no_dbWB)=frequency(iif);
                glm_dB_hipp.per_ii(glm_ii_hipp+1:glm_ii_hipp+no_dbWB)=per_ii;
                glm_dB_hipp.evNo(glm_ii_hipp+1:glm_ii_hipp+no_dbWB)=evNo;
                glm_dB_hipp.grNo(glm_ii_hipp+1:glm_ii_hipp+no_dbWB)=grNo;
                glm_ii_hipp=glm_ii_hipp+no_dbWB;
                
                glm_dB_both.data(glm_ii_both+1:glm_ii_both+no_dbWB)=these_dBWB(:,iif);
                glm_dB_both.freq(glm_ii_both+1:glm_ii_both+no_dbWB)=frequency(iif);
                glm_dB_both.per_ii(glm_ii_both+1:glm_ii_both+no_dbWB)=per_ii;
                glm_dB_both.evNo(glm_ii_both+1:glm_ii_both+no_dbWB)=evNo;
                glm_dB_both.grNo(glm_ii_both+1:glm_ii_both+no_dbWB)=grNo;
                glm_dB_both.tissue(glm_ii_both+1:glm_ii_both+no_dbWB)=1;
                glm_ii_both=glm_ii_both+no_dbWB;
            end
            
            maxdB=max([maxdB max(this_mean_dbWB)]);
            mindB=min([mindB min(this_mean_dbWB)]);
            
            
        end
    end
    
    title(['Delta LFP power vs frequency for hippocampus' group_legend{grNo}])
    
    
    %Proficient/Naive annotations
    annotation('textbox',[0.15 0.85 0.3 0.1],'String','S+ Proficient','FitBoxToText','on','Color',[158/255 31/255 99/255],'LineStyle','none');
    annotation('textbox',[0.15 0.80 0.3 0.1],'String','S+ Naive','FitBoxToText','on','Color',[238/255 111/255 179/255],'LineStyle','none');
    annotation('textbox',[0.15 0.75 0.3 0.1],'String','S- Proficient','FitBoxToText','on','Color',[[0 114/255 178/255]],'LineStyle','none');
    annotation('textbox',[0.15 0.70 0.3 0.1],'String','S- Naive','FitBoxToText','on','Color',[80/255 194/255 255/255],'LineStyle','none');
    
    
    xlim([4 95])
    
    xlabel('Frequency (Hz)')
    
    ylim([mindB-0.05*(maxdB-mindB) maxdB+0.05*(maxdB-mindB)])
    
    
    ylabel('dB')
    
    
    
    %Now plot for the prefrontal the dependence of delta LFP to frequency
    
    %Plot the average
    figNo = figNo +1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    
    
    set(hFig, 'units','normalized','position',[.1 .5 .35 .3])
    hold on
    
    
    for evNo=1:2
        
        for per_ii=2:-1:1
            
            
            
            %Get these MI values
            these_dBWB=[];
            ii_dB=0;
            for ii=1:length(hippFileName)
                this_jj=[];
                for jj=1:length(all_pre(ii).handles_out.delta_dB_WB)
                    if all_pre(ii).handles_out.delta_dB_WB(jj).per_ii==per_ii
                        if all_pre(ii).handles_out.delta_dB_WB(jj).evNo==evNo
                            if all_pre(ii).handles_out.delta_dB_WB(jj).grNo==grNo
                                this_jj=jj;
                            end
                        end
                    end
                end
                if ~isempty(this_jj)
                    ii_dB=ii_dB+1;
                    these_dBWB(ii_dB,1:length(frequency))=all_pre(ii).handles_out.delta_dB_WB(this_jj).mean_dbWB;
                end
            end
            
            this_mean_dbWB=zeros(1,length(frequency));
            this_mean_dbWB(1,:)=mean(these_dBWB,1);
            
            CI=[];
            CI = bootci(1000, {@mean, these_dBWB})';
            CI(:,1)= this_mean_dbWB'-CI(:,1);
            CI(:,2)=CI(:,2)- this_mean_dbWB';
            
            if evNo==1
                if per_ii==1
                    %S+ Proficient
                    [hlCR, hpCR] = boundedline(frequency',this_mean_dbWB', CI, 'cmap',[158/255 31/255 99/255]);
                else
                    %S+ Naive
                    [hlCR, hpCR] = boundedline(frequency',this_mean_dbWB', CI, 'cmap',[238/255 111/255 179/255]);
                end
            else
                if per_ii==1
                    %S- Proficient
                    [hlCR, hpCR] = boundedline(frequency',this_mean_dbWB', CI, 'cmap',[0 114/255 178/255]);
                else
                    %S- naive
                    [hlCR, hpCR] = boundedline(frequency',this_mean_dbWB', CI, 'cmap',[80/255 194/255 255/255]);
                end
            end
            
            no_dbWB=size(these_dBWB,1);
            for iif=1:length(frequency)
                glm_dB_pre.data(glm_ii_pre+1:glm_ii_pre+no_dbWB)=these_dBWB(:,iif);
                glm_dB_pre.freq(glm_ii_pre+1:glm_ii_pre+no_dbWB)=frequency(iif);
                glm_dB_pre.per_ii(glm_ii_pre+1:glm_ii_pre+no_dbWB)=per_ii;
                glm_dB_pre.evNo(glm_ii_pre+1:glm_ii_pre+no_dbWB)=evNo;
                glm_dB_pre.grNo(glm_ii_pre+1:glm_ii_pre+no_dbWB)=grNo;
                glm_ii_pre=glm_ii_pre+no_dbWB;
                
                glm_dB_both.data(glm_ii_both+1:glm_ii_both+no_dbWB)=these_dBWB(:,iif);
                glm_dB_both.freq(glm_ii_both+1:glm_ii_both+no_dbWB)=frequency(iif);
                glm_dB_both.per_ii(glm_ii_both+1:glm_ii_both+no_dbWB)=per_ii;
                glm_dB_both.grNo(glm_ii_both+1:glm_ii_both+no_dbWB)=grNo;
                glm_dB_both.evNo(glm_ii_both+1:glm_ii_both+no_dbWB)=evNo;
                glm_dB_both.tissue(glm_ii_both+1:glm_ii_both+no_dbWB)=2;
                glm_ii_both=glm_ii_both+no_dbWB;
            end
            
            maxdB=max([maxdB max(this_mean_dbWB)]);
            mindB=min([mindB min(this_mean_dbWB)]);
            
            
        end
    end
    
    title(['Delta LFP power vs frequency for prefrontal' group_legend{grNo}])
    
    
    %Proficient/Naive annotations
    annotation('textbox',[0.15 0.85 0.3 0.1],'String','S+ Proficient','FitBoxToText','on','Color',[158/255 31/255 99/255],'LineStyle','none');
    annotation('textbox',[0.15 0.80 0.3 0.1],'String','S+ Naive','FitBoxToText','on','Color',[238/255 111/255 179/255],'LineStyle','none');
    annotation('textbox',[0.15 0.75 0.3 0.1],'String','S- Proficient','FitBoxToText','on','Color',[[0 114/255 178/255]],'LineStyle','none');
    annotation('textbox',[0.15 0.70 0.3 0.1],'String','S- Naive','FitBoxToText','on','Color',[80/255 194/255 255/255],'LineStyle','none');
    
    
    xlim([4 95])
    
    xlabel('Frequency (Hz)')
    
    ylim([mindB-0.05*(maxdB-mindB) maxdB+0.05*(maxdB-mindB)])
    
    
    ylabel('dB')
    
end

for fn=1:figNo
    figure(fn)
    ylim([mindB-0.05*(maxdB-mindB) maxdB+0.05*(maxdB-mindB)])
end

%Perform the glm
fprintf(1, ['\n\nglm for delta dB calculated per mouse for hippocampus\n'])
tbl = table(glm_dB_hipp.data',glm_dB_hipp.freq',glm_dB_hipp.evNo',glm_dB_hipp.per_ii',glm_dB_hipp.grNo',...
    'VariableNames',{'dB','frequency','spm','naive_proficient','group'});
mdl = fitglm(tbl,'dB~frequency+spm+naive_proficient+group+frequency*spm*naive_proficient*group'...
    ,'CategoricalVars',[3,4,5])
 
%Perform the glm
fprintf(1, ['\n\nglm for delta dB calculated per mouse for prefrontal\n'])
tbl = table(glm_dB_pre.data',glm_dB_pre.freq',glm_dB_pre.evNo',glm_dB_pre.per_ii',glm_dB_pre.grNo',...
    'VariableNames',{'dB','frequency','spm','naive_proficient','group'});
mdl = fitglm(tbl,'dB~frequency+spm+naive_proficient+group+frequency*spm*naive_proficient*group'...
    ,'CategoricalVars',[3,4,5])


%Perform the glm mi for both brain regions
fprintf(1, ['\n\nglm for delta dB for prefrontal and hippocampus\n'])
tbl = table(glm_dB_both.data',glm_dB_both.freq',glm_dB_both.evNo',glm_dB_both.per_ii',glm_dB_both.grNo',glm_dB_both.tissue',...
    'VariableNames',{'dB','frequency','spm','naive_proficient','group','tissue'});
mdl = fitglm(tbl,'dB~frequency+spm+naive_proficient+group+tissue+frequency*spm*naive_proficient*group*tissue'...
    ,'CategoricalVars',[3,4,5,6])



pffft=1

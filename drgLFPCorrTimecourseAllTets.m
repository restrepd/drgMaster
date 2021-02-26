function handles=drgLFPCorrTimecourseAllTets(handles)

handles.displayData=0;
handles.save_drgb=1;
handles_corr_out=[];

if ispc
    handles.tempPath=[handles.jtPathName '\temp_corr\'];
else
    handles.tempPath=[handles.jtPathName '/temp_corr/'];
end

if ~isfield(handles,'recalculate')
    handles.recalc_corr=0;
end

if ~isfield(handles,'drgbchoices')
    handles.drgbchoices.other_electrodes=[9:16]; %Prefrontal
    handles.drgbchoices.other_legend='Prefrontal';
    handles.drgbchoices.reference_electrodes=[1:8]; %Hippocampus
    handles.drgbchoices.reference_legend='Hippocampus';
    %5 is S+ and 11 is S-
    handles.drgbchoices.evTypeNos=[3 5 7 9 11 13]; 
    handles.drgbchoices.referenceEvent=2;
end

if (exist([handles.jtPathName 'Corr_' handles.jtFileName(9:end)])~=2)||...
        (handles.recalc_corr==1)
    gcp;
    
    %handles.burstLFPNo reference LFP (hippocampus)
    %handles.peakLFPNo the other LFP
    
    %Between tetrodes
    handles_corr_out.no_between_pairs=0;
    par_out_LFPcorr=[];
    
    
    
    for ref_tet_electNo=handles.drgbchoices.reference_electrodes
        handles.burstLFPNo=ref_tet_electNo;
        for other_tet_electNo=handles.drgbchoices.other_electrodes
            handles.peakLFPNo=other_tet_electNo;
            handles_corr_out.no_between_pairs=handles_corr_out.no_between_pairs+1;
            handles_corr_out.between.burstLFPNo(handles_corr_out.no_between_pairs)=handles.burstLFPNo;
            handles_corr_out.between.peakLFPNo(handles_corr_out.no_between_pairs)=handles.peakLFPNo;
            par_out_LFPcorr(handles_corr_out.no_between_pairs).LFPcorr=[];
        end
    end
    
    
    %Run drgLFPCorrTimecourse
    no_between_pairs=handles_corr_out.no_between_pairs;
    parfor ii_between_pairs=1:no_between_pairs
        % for ii_between_pairs=1:no_between_pairs
        handlespf=struct();
        handlespf=handles;
        handlespf.corr_out=handles_corr_out;
        handlespf.burstLFPNo=handlespf.corr_out.between.burstLFPNo(ii_between_pairs);
        handlespf.peakLFPNo=handlespf.corr_out.between.peakLFPNo(ii_between_pairs);
        
        %Do between tetrode calculation
        handlespf.randpermLFP=0;
        if handlespf.recalc_corr==1
            handlespf=drgLFPCorrTimecourse(handlespf);
            par_out_LFPcorr(ii_between_pairs).LFPcorr=handlespf.drgb.LFPcorr;
            this_LFPcorr=handlespf.drgb.LFPcorr;
            drgSaveParCorr(handlespf.tempPath,['temp_corr_between' num2str(ii_between_pairs) '.mat'],this_LFPcorr)
        else
            if exist([handlespf.tempPath 'temp_corr_between' num2str(ii_between_pairs) '.mat'],'file')==2
                par_out_LFPcorr(ii_between_pairs).LFPcorr=drgLoadParCorr(handlespf.tempPath,['temp_corr_between' num2str(ii_between_pairs) '.mat']);
            else
                handlespf=drgLFPCorrTimecourse(handlespf);
                par_out_LFPcorr(ii_between_pairs).LFPcorr=handlespf.drgb.LFPcorr;
                this_LFPcorr=handlespf.drgb.LFPcorr;
                drgSaveParCorr(handlespf.tempPath,['temp_corr_between' num2str(ii_between_pairs) '.mat'],this_LFPcorr)
            end
        end
        fprintf(1, ['Between brain regions pair No %d, reference electrode %d, other electrode %d\n'], ii_between_pairs, handlespf.burstLFPNo,handlespf.peakLFPNo);
        
%         %Do between tetrode control calculation by flipping one LFP
%         handlespf.randpermLFP=1;
%         if handlespf.recalc_corr==1
%             handlespf=drgLFPCorrTimecourse(handlespf);
%             par_out_LFPcorr(ii_between_pairs).LFPcorr_randperm=handlespf.drgb.LFPcorr;
%             this_LFPcorr=handlespf.drgb.LFPcorr;
%             drgSaveParCorr(handlespf.tempPath,['temp_corr_randperm_between' num2str(ii_between_pairs) '.mat'],this_LFPcorr)
%         else
%             if exist([handlespf.tempPath 'temp_corr_randperm_between' num2str(ii_between_pairs) '.mat'],'file')==2
%                 par_out_LFPcorr(ii_between_pairs).LFPcorr_randperm=drgLoadParCorr(handlespf.tempPath,['temp_corr_randperm_between' num2str(ii_between_pairs) '.mat']);
%             else
%                 handlespf=drgLFPCorrTimecourse(handlespf);
%                 par_out_LFPcorr(ii_between_pairs).LFPcorr_randperm=handlespf.drgb.LFPcorr;
%                 this_LFPcorr=handlespf.drgb.LFPcorr;
%                 drgSaveParCorr(handlespf.tempPath,['temp_corr_randperm_between' num2str(ii_between_pairs) '.mat'],this_LFPcorr)
%             end
%         end
%         fprintf(1, ['Between brain regions pair No %d, reference electrode %d, other electrode, %d LFP permutation\n'], ii_between_pairs, handlespf.burstLFPNo,handlespf.peakLFPNo);
    end
    
    %Save the output to handles
    for ii_between_pairs=1:no_between_pairs
        handles_corr_out.between.LFPcorr(ii_between_pairs).LFPcorr=par_out_LFPcorr(ii_between_pairs).LFPcorr;
%         handles_corr_out.between.LFPcorr_randperm(ii_between_pairs).LFPcorr=par_out_LFPcorr(ii_between_pairs).LFPcorr_randperm;
    end
    
    
    %Within tetrodes
    handles_corr_out.no_within_pairs=0;
    par_out_LFPcorr=[];
    
    
    
    for ref_tet_electNo1=handles.drgbchoices.reference_electrodes
        handles.burstLFPNo=ref_tet_electNo1;
        this_el=find(handles.drgbchoices.reference_electrodes==ref_tet_electNo1);
        for ref_tet_electNo2=handles.drgbchoices.reference_electrodes(this_el+1:end)
            handles.peakLFPNo=ref_tet_electNo2;
            handles_corr_out.no_within_pairs=handles_corr_out.no_within_pairs+1;
            handles_corr_out.within.burstLFPNo(handles_corr_out.no_within_pairs)=handles.burstLFPNo;
            handles_corr_out.within.peakLFPNo(handles_corr_out.no_within_pairs)=handles.peakLFPNo;
%             handles_corr_out.within.refTetNo(handles_corr_out.no_within_pairs)=ref_tetNo;
            par_out_LFPcorr(handles_corr_out.no_within_pairs).LFPcorr=[];
        end
    end
    
    for ref_tet_electNo1=handles.drgbchoices.other_electrodes
        handles.burstLFPNo=ref_tet_electNo1;
        this_el=find(handles.drgbchoices.other_electrodes==ref_tet_electNo1);
        for ref_tet_electNo2=handles.drgbchoices.other_electrodes(this_el+1:end)
            handles.peakLFPNo=ref_tet_electNo2;
            handles_corr_out.no_within_pairs=handles_corr_out.no_within_pairs+1;
            handles_corr_out.within.burstLFPNo(handles_corr_out.no_within_pairs)=handles.burstLFPNo;
            handles_corr_out.within.peakLFPNo(handles_corr_out.no_within_pairs)=handles.peakLFPNo;
%             handles_corr_out.within.refTetNo(handles_corr_out.no_within_pairs)=ref_tetNo;
            par_out_LFPcorr(handles_corr_out.no_within_pairs).LFPcorr=[];
        end
    end
    
    
    %Run drgLFPCorrTimecourse
    no_within_pairs=handles_corr_out.no_within_pairs;
    parfor ii_within_pairs=1:no_within_pairs
        % for ii_within_pairs=1:no_within_pairs
        handlespf=struct();
        handlespf=handles;
        handlespf.corr_out=handles_corr_out;
        handlespf.burstLFPNo=handlespf.corr_out.within.burstLFPNo(ii_within_pairs);
        handlespf.peakLFPNo=handlespf.corr_out.within.peakLFPNo(ii_within_pairs);
        
        %Do within tetrode calculation
        handlespf.randpermLFP=0;
        if handlespf.recalc_corr==1
            handlespf=drgLFPCorrTimecourse(handlespf);
            par_out_LFPcorr(ii_within_pairs).LFPcorr=handlespf.drgb.LFPcorr;
            this_LFPcorr=handlespf.drgb.LFPcorr;
            drgSaveParCorr(handlespf.tempPath,['temp_corr_within' num2str(ii_within_pairs) '.mat'],this_LFPcorr)
        else
            if exist([handlespf.tempPath 'temp_corr_within' num2str(ii_within_pairs) '.mat'],'file')==2
                par_out_LFPcorr(ii_within_pairs).LFPcorr=drgLoadParCorr(handlespf.tempPath,['temp_corr_within' num2str(ii_within_pairs) '.mat']);
            else
                handlespf=drgLFPCorrTimecourse(handlespf);
                par_out_LFPcorr(ii_within_pairs).LFPcorr=handlespf.drgb.LFPcorr;
                this_LFPcorr=handlespf.drgb.LFPcorr;
                drgSaveParCorr(handlespf.tempPath,['temp_corr_within' num2str(ii_within_pairs) '.mat'],this_LFPcorr)
            end
        end
        fprintf(1, ['Within brain region pair No %d, reference electrode %d, other electrode %d\n'], ii_within_pairs, handlespf.burstLFPNo,handlespf.peakLFPNo);
        
%         %Do within tetrode control calculation by flipping one LFP
%         handlespf.randpermLFP=1;
%         if handlespf.recalc_corr==1
%             handlespf=drgLFPCorrTimecourse(handlespf);
%             par_out_LFPcorr(ii_within_pairs).LFPcorr_randperm=handlespf.drgb.LFPcorr;
%             this_LFPcorr=handlespf.drgb.LFPcorr;
%             drgSaveParCorr(handlespf.tempPath,['temp_corr_randperm_within' num2str(ii_within_pairs) '.mat'],this_LFPcorr)
%         else
%             if exist([handlespf.tempPath 'temp_corr_randperm_within' num2str(ii_within_pairs) '.mat'],'file')==2
%                 par_out_LFPcorr(ii_within_pairs).LFPcorr=drgLoadParCorr(handlespf.tempPath,['temp_corr_within' num2str(ii_within_pairs) '.mat']);
%             else
%                 handlespf=drgLFPCorrTimecourse(handlespf);
%                 par_out_LFPcorr(ii_within_pairs).LFPcorr_randperm=handlespf.drgb.LFPcorr;
%                 this_LFPcorr=handlespf.drgb.LFPcorr;
%                 drgSaveParCorr(handlespf.tempPath,['temp_corr_randperm_within' num2str(ii_within_pairs) '.mat'],this_LFPcorr)
%             end
%         end
%         fprintf(1, ['Within tetrodes pair No %d, reference electrode %d, other electrode %d, LFP permutation\n'], ii_within_pairs, handlespf.burstLFPNo,handlespf.peakLFPNo);
    end
    
    %Save the output to handles
    for ii_within_pairs=1:no_within_pairs
        handles_corr_out.within.LFPcorr(ii_within_pairs).LFPcorr=par_out_LFPcorr(ii_within_pairs).LFPcorr;
%         handles_corr_out.within.LFPcorr_randperm(ii_within_pairs).LFPcorr=par_out_LFPcorr(ii_within_pairs).LFPcorr_randperm;
    end
    
    %Turn off the parallel pool and release RAM
    poolobj = gcp('nocreate');
    delete(poolobj);
    
    %Save the cross-correlation data
    save([handles.jtPathName 'Corr_' handles.jtFileName(9:end)],'handles_corr_out','-v7.3')
    
else
    if (exist([handles.jtPathName 'Corr_' handles.jtFileName(9:end)])==2)
        load([handles.jtPathName 'Corr_' handles.jtFileName(9:end)],'handles_corr_out')
    end
end

%Now plot the max_rho_t_lag
try
    close 3
catch
end


hFig3 = figure(3);
set(hFig3, 'units','normalized','position',[.26 .05 .3 .3])

hold on



no_within_pairs=length(handles_corr_out.within.LFPcorr);
no_between_pairs=length(handles_corr_out.between.LFPcorr);

% 
% % subplot(2,1,2)
% % hold on
% 
% all_max_rho_t_lag=[];
% all_mean_max_rho_t_lag=[];
% ii=0;
% time=handles_corr_out.within.LFPcorr_randperm(1).LFPcorr.time;
% for ii_within_pairs=1:no_within_pairs
%     
%     max_rho_t_lag=zeros(length(handles_corr_out.within.LFPcorr_randperm(ii_within_pairs).LFPcorr.trial),length(time));
%     for trNo=1:length(handles_corr_out.within.LFPcorr_randperm(ii_within_pairs).LFPcorr.trial)
%         max_rho_t_lag(trNo,:)=handles_corr_out.within.LFPcorr_randperm(ii_within_pairs).LFPcorr.trial(trNo).max_rho_t_lag;
%     end
%     mean_max_rho_t_lag=mean(max_rho_t_lag);
%     ii=ii+1;
%     all_max_rho_t_lag=[all_max_rho_t_lag; max_rho_t_lag];
%     all_mean_max_rho_t_lag(ii,:)=mean_max_rho_t_lag;
%     
% end
% 
% 
% 
% mean_all_max_rho_t_lag=mean(all_mean_max_rho_t_lag,1);
% tempCI = bootci(1000, {@mean, all_mean_max_rho_t_lag})';
% CI=[];
% CI(:,1)=mean_all_max_rho_t_lag'-tempCI(:,1);
% CI(:,2)=tempCI(:,2)-mean_all_max_rho_t_lag';
% 
% [hlCR, hpCR] = boundedline(time,mean_all_max_rho_t_lag, CI, 'b');

spm=[2 5];
pcorr=[45 65;80 100];
spm_legends={'S+','S-'};
pcorr_legends={'Naive','Proficient'};
spii=0;



for ii_pcorr=1:2
    for ii_spm=1:2
        spii=spii+1;
        subplot(2,2,spii)
        hold on
        
%         %Plot within pair
%         all_max_rho_t_lag=[];
%         all_mean_max_rho_t_lag=[];
%         ii=0;
%         time=handles_corr_out.within.LFPcorr(1).LFPcorr.time;
%         for ii_within_pairs=1:no_within_pairs
%             these_trials=(handles_corr_out.within.LFPcorr(ii_within_pairs).LFPcorr.which_event(spm(ii_spm),:))&...
%                 (handles_corr_out.within.LFPcorr(ii_within_pairs).LFPcorr.perCorr>=pcorr(ii_pcorr,1))&...
%                 (handles_corr_out.within.LFPcorr(ii_within_pairs).LFPcorr.perCorr<=pcorr(ii_pcorr,2));
%             max_rho_t_lag=zeros(sum(these_trials),length(time));
%             no_trials=0;
%             for trNo=1:length(handles_corr_out.within.LFPcorr(ii_within_pairs).LFPcorr.trial)
%                 if these_trials(trNo)==1
%                     no_trials=no_trials+1;
%                     max_rho_t_lag(no_trials,:)=handles_corr_out.within.LFPcorr(ii_within_pairs).LFPcorr.trial(trNo).max_rho_t_lag;
%                 end
%             end
%             mean_max_rho_t_lag=mean(max_rho_t_lag);
%             ii=ii+1;
%             all_max_rho_t_lag=[all_max_rho_t_lag; max_rho_t_lag];
%             all_mean_max_rho_t_lag(ii,:)=mean_max_rho_t_lag;
%         end
%         
%         
%         
%         mean_all_max_rho_t_lag=mean(all_mean_max_rho_t_lag,1);
%         tempCI = bootci(1000, {@mean, all_mean_max_rho_t_lag})';
%         CI=[];
%         CI(:,1)=mean_all_max_rho_t_lag'-tempCI(:,1);
%         CI(:,2)=tempCI(:,2)-mean_all_max_rho_t_lag';
%         
%         [hlCR, hpCR] = boundedline(time,mean_all_max_rho_t_lag, CI, 'b');
        
        %Plot between pair graph
        
        
        
        % subplot(2,1,1)
        % hold on
        %
        % all_max_rho_t_lag=[];
        % all_mean_max_rho_t_lag=[];
        % ii=0;
        % time=handles_corr_out.between.LFPcorr_randperm(1).LFPcorr.time;
        % for ii_between_pairs=1:no_between_pairs
        %
        %     max_rho_t_lag=zeros(length(handles_corr_out.between.LFPcorr_randperm(ii_between_pairs).LFPcorr.trial),length(time));
        %     for trNo=1:length(handles_corr_out.between.LFPcorr_randperm(ii_between_pairs).LFPcorr.trial)
        %         max_rho_t_lag(trNo,:)=handles_corr_out.between.LFPcorr_randperm(ii_between_pairs).LFPcorr.trial(trNo).max_rho_t_lag;
        %     end
        %     mean_max_rho_t_lag=mean(max_rho_t_lag);
        %     ii=ii+1;
        %     all_max_rho_t_lag=[all_max_rho_t_lag; max_rho_t_lag];
        %     all_mean_max_rho_t_lag(ii,:)=mean_max_rho_t_lag;
        %
        % end
        %
        %
        %
        % mean_all_max_rho_t_lag=mean(all_mean_max_rho_t_lag,1);
        % tempCI = bootci(1000, {@mean, all_mean_max_rho_t_lag})';
        % CI=[];
        % CI(:,1)=mean_all_max_rho_t_lag'-tempCI(:,1);
        % CI(:,2)=tempCI(:,2)-mean_all_max_rho_t_lag';
        %
        % [hlCR, hpCR] = boundedline(time,mean_all_max_rho_t_lag, CI, 'b');
        
        all_max_rho_t_lag=[];
        all_mean_max_rho_t_lag=[];
        ii=0;
        time=handles_corr_out.between.LFPcorr(1).LFPcorr.time;
        for ii_between_pairs=1:no_between_pairs
            these_trials=(handles_corr_out.between.LFPcorr(ii_between_pairs).LFPcorr.which_event(spm(ii_spm),:))&...
                (handles_corr_out.between.LFPcorr(ii_between_pairs).LFPcorr.perCorr>=pcorr(ii_pcorr,1))&...
                (handles_corr_out.between.LFPcorr(ii_between_pairs).LFPcorr.perCorr<=pcorr(ii_pcorr,2));
            max_rho_t_lag=zeros(sum(these_trials),length(time));
            no_trials=0;
            for trNo=1:length(handles_corr_out.between.LFPcorr(ii_between_pairs).LFPcorr.trial)
                if these_trials(trNo)==1
                    no_trials=no_trials+1;
                    max_rho_t_lag(no_trials,:)=handles_corr_out.between.LFPcorr(ii_between_pairs).LFPcorr.trial(trNo).max_rho_t_lag;
                end
            end
            
            %             max_rho_t_lag=zeros(length(handles_corr_out.between.LFPcorr(ii_between_pairs).LFPcorr.trial),length(time));
            %             for trNo=1:length(handles_corr_out.between.LFPcorr(ii_between_pairs).LFPcorr.trial)
            %                 max_rho_t_lag(trNo,:)=handles_corr_out.between.LFPcorr(ii_between_pairs).LFPcorr.trial(trNo).max_rho_t_lag;
            %             end
            mean_max_rho_t_lag=mean(max_rho_t_lag);
            ii=ii+1;
            all_max_rho_t_lag=[all_max_rho_t_lag; max_rho_t_lag];
            all_mean_max_rho_t_lag(ii,:)=mean_max_rho_t_lag;
            
        end
        
        
        
        mean_all_max_rho_t_lag=mean(all_mean_max_rho_t_lag,1);
        tempCI = bootci(1000, {@mean, all_mean_max_rho_t_lag})';
        CI=[];
        CI(:,1)=mean_all_max_rho_t_lag'-tempCI(:,1);
        CI(:,2)=tempCI(:,2)-mean_all_max_rho_t_lag';
        
        [hlCR, hpCR] = boundedline(time,mean_all_max_rho_t_lag, CI, 'r');
        
        hold on
        plot([time(1) time(end)],[0 0])
        ylim([-0.02 0.02])
        this_ylm=ylim;
        title([spm_legends{ii_spm} ' ' pcorr_legends{ii_pcorr}])
        %         text(-0.8,0.9*(this_ylm(2)-this_ylm(1))+this_ylm(1),['Between, reference  ' handles.drgbchoices.reference_legend ' vs ' handles.drgbchoices.other_legend],'Color','r')
        %         text(-0.8,0.85*(this_ylm(2)-this_ylm(1))+this_ylm(1),'Within ','Color','b')
        xlabel('Time (sec)')
        ylabel('Lag (sec)')
        
    end
end

suptitle(['Max rho lag timecourse reference  ' handles.drgbchoices.reference_legend ' vs ' handles.drgbchoices.other_legend])

%Now plot a pseudocolor for the mean lag
try
    close 1
catch
end

hFig1 = figure(1);
set(hFig1, 'units','normalized','position',[.05 .05 .2 .8])

max_rho_t_lag=zeros(sum(these_trials),length(time),no_between_pairs);

for ii_between_pairs=1:no_between_pairs
    for trNo=1:length(handles_corr_out.between.LFPcorr(ii_between_pairs).LFPcorr.trial) 
        max_rho_t_lag(trNo,:,ii_between_pairs)=handles_corr_out.between.LFPcorr(ii_between_pairs).LFPcorr.trial(trNo).max_rho_t_lag;   
    end
end

mean_max_rho_t_lag=mean(max_rho_t_lag,3);

min_prob=prctile(mean_max_rho_t_lag(:),5);
max_prob=prctile(mean_max_rho_t_lag(:),95);


%pcolor(repmat(phase,length(trials),1),repmat(trials',1,length(phase)),all_phase_histo)
trials=1:length(handles_corr_out.between.LFPcorr(ii_between_pairs).LFPcorr.trial);
no_trials=length(handles_corr_out.between.LFPcorr(ii_between_pairs).LFPcorr.trial);
drg_pcolor(repmat(time,no_trials,1),repmat(trials',1,length(time)),mean_max_rho_t_lag)
%         colormap jet
colormap fire
shading flat
% min_prob=0.0113;
% max_prob=0.0314;
caxis([min_prob    max_prob])
xlabel('Time (sec)')
ylabel('Trial');
title(['Lag (sec)'])


try
    close 2
catch
end

hFig2 = figure(2);
set(hFig2, 'units','normalized','position',[.26 .5 .05 .3])


prain=[min_prob:(max_prob-min_prob)/99:max_prob];
pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
%             colormap jet
colormap fire
shading interp
ax=gca;
set(ax,'XTickLabel','')


pffft=1;
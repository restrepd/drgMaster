close all
clear all

[PLVFileName,PLVBatchPathName] = uigetfile({'*.mat'},'Select the .mat file with the results of drgRunBatchPLVpar');
fprintf(1, ['\ndrgAnalyzePLVBatcg run for ' PLVFileName '\n\n']);
 
load([PLVBatchPathName PLVFileName])

%Display the data
[choiceFileName,choiceBatchPathName] = uigetfile({'drgbChoicesPLV*.m'},'Select the .m file with all the choices for analysis');
        
tempDirName=['temp' choiceFileName(12:end-2)];

addpath(choiceBatchPathName)
eval(['handles_in=' choiceFileName(1:end-2) ';'])

handles.drgbchoices=handles_in.drgbchoices;

figNo=0;
 


for grNo=1:3
           
    
    %Plot average PLV timecourse between brain regions
            
            
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            
            
            hFig = figure(figNo);
            set(hFig, 'units','normalized','position',[.05 .05 .45 .7])
            
            this_sub=0;
            
            for bw_ii=1:length(handles.drgbchoices.lowF)
                
                for per_ii=size(handles.drgbchoices.percent_windows,1):-1:1
                    this_sub=this_sub+1;
                    subplot(4,2,this_sub)
                     subplot(length(handles.drgbchoices.lowF),size(handles.drgbchoices.percent_windows,1),this_sub)
                    hold on
                    for eventNo_ii=length(handles.drgbchoices.evTypeNos):-1:1
                        
                        %Get the PLVs
                        these_ii_PLVs=[];
                        jj=0;
                        
                        for ii=1:length(handles_out.PLV)
                            if (handles_out.PLV(ii).eventNo_ii==eventNo_ii)&(handles_out.PLV(ii).per_ii==per_ii)&(handles_out.PLV(ii).bw_ii==bw_ii)&(handles_out.PLV(ii).group_no==grNo)
                                jj=jj+1;
                                these_ii_PLVs(jj)=ii;
                            end
                        end
                        
                        mean_plvs=[];
                        if ~isfield(handles_out,'drg')
                            Fs=20000;
                        else
                            Fs=handles_out.drg.session(1).draq_p.ActualRate;
                        end
                        
                        for ii_plv_per_mouse=1:length(these_ii_PLVs)
                            this_ii_PLV=these_ii_PLVs(ii_plv_per_mouse);
                            these_PLVs=zeros(length(handles.drgbchoices.reference_electrodes)*length(handles.drgbchoices.other_electrodes),length(handles_out.PLV(this_ii_PLV).time));
                            ii_PLVs=0;
                            for ii_ref=handles.drgbchoices.reference_electrodes
                                for ii_oth=handles.drgbchoices.other_electrodes
                                    this_plv=zeros(1,length(handles_out.PLV(this_ii_PLV).time));
                                    if ii_ref<ii_oth
                                        this_plv(1,:)=handles_out.PLV(this_ii_PLV).plv(:,ii_ref,ii_oth);
                                    else
                                        this_plv(1,:)=handles_out.PLV(this_ii_PLV).plv(:,ii_oth,ii_ref);
                                    end
                                    ii_PLVs=ii_PLVs+1;
                                    these_PLVs(ii_PLVs,:)=this_plv;
                                end
                            end
                            
                            %Calculate plv for 100 msec intervals
                            dt=0.1;
                            decimation_factor=40;
                            no_time_points=(Fs/decimation_factor)*dt;
                            these_PLVs_dec=zeros(size(these_PLVs,1),size(these_PLVs,2)/no_time_points);
                            time_dec=zeros(1,size(these_PLVs,2)/no_time_points);
                            time=handles_out.PLV(this_ii_PLV).time;
                            
                            for elec_pair=1:size(these_PLVs,1)
                                for ii_dec_tp=1:size(these_PLVs,2)/no_time_points
                                    these_PLVs_dec(elec_pair,ii_dec_tp)=mean(these_PLVs(elec_pair,1+(ii_dec_tp-1)*no_time_points:ii_dec_tp*no_time_points));
                                    if elec_pair==1
                                        time_dec(1,ii_dec_tp)=mean(time(1+(ii_dec_tp-1)*no_time_points:ii_dec_tp*no_time_points));
                                    end
                                end
                            end
                            
                            
                            mean_plv(ii_plv_per_mouse,:)=mean(these_PLVs_dec);
                            
                        end
                        
%                         CI = bootci(1000, {@mean, mean_plv})';
%                         CI(:,1)=mean(mean_plv)'-CI(:,1);
%                         CI(:,2)=CI(:,2)-mean(mean_plv)';
%                         
%                         
%                         if eventNo_ii==2
%                             [hlCR, hpCR] = boundedline(time_dec,mean(mean_plv), CI, 'b');
%                         else
%                             [hlCR, hpCR] = boundedline(time_dec,mean(mean_plv), CI, 'r');
%                         end

           
                        if eventNo_ii==2
                            plot(time_dec,mean(mean_plv), 'b');
                        else
                            plot(time_dec,mean(mean_plv), 'r');
                        end
                        
                        
                        
                        %Get PLS
                        these_PLVs=zeros(length(handles.drgbchoices.reference_electrodes),length(handles_out.PLV(this_ii_PLV).time));
                        ii_PLVs=0;
                        for ii_ref=1:handles.drgbchoices.no_electrodes
                            try
                                ii_oth=handles.drgbchoices.no_electrodes+1;
                                this_plv=zeros(1,length(handles_out.PLV(this_ii_PLV).time));
                                if ii_ref<ii_oth
                                    this_plv(1,:)=handles_out.PLV(this_ii_PLV).plv(:,ii_ref,ii_oth);
                                else
                                    this_plv(1,:)=handles_out.PLV(this_ii_PLV).plv(:,ii_oth,ii_ref);
                                end
                            catch
                            end
                            ii_PLVs=ii_PLVs+1;
                            these_PLVs(ii_PLVs,:)=this_plv;
                            
                        end
                        pls=prctile(these_PLVs(:),95);
                        plot([time_dec(1) time_dec(end)],[pls pls],'-k')
                        
                        ylim([0 1])
                        xlabel('Time (sec)')
                        ylabel('PLV')
                        title([handles.drgbchoices.bwlabels{bw_ii} ' ' handles.drgbchoices.percent_labels{per_ii}])
                        
                        if (bw_ii==1)&(per_ii==size(handles.drgbchoices.percent_windows,1))
                            text(0.95,-1.5,'S+','Color','r')
                            text(0.9,-1.5,'S-','Color','b')
                            text(0.9,-1.5,'p<0.05','Color','k')
                        end
                        
                    end
                end
            end
            
            suptitle(['Phase-locking value for ' handles.drgbchoices.group_no_names{grNo}])
            
            
            
            %Plot delta phase timecourse between brain regions
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            
            
            hFig = figure(figNo);
            set(hFig, 'units','normalized','position',[.05 .05 .45 .7])
            
            this_sub=0;
            
            for bw_ii=1:length(handles.drgbchoices.lowF)
                
                for per_ii=size(handles.drgbchoices.percent_windows,1):-1:1
                    this_sub=this_sub+1;
                    subplot(length(handles.drgbchoices.lowF),size(handles.drgbchoices.percent_windows,1),this_sub)
                    hold on
                    for eventNo_ii=length(handles.drgbchoices.evTypeNos):-1:1
                        
                       %Get the PLVs
                        these_ii_PLVs=[];
                        jj=0;
                        
                        for ii=1:length(handles_out.PLV)
                            if (handles_out.PLV(ii).eventNo_ii==eventNo_ii)&(handles_out.PLV(ii).per_ii==per_ii)&(handles_out.PLV(ii).bw_ii==bw_ii)&(handles_out.PLV(ii).group_no==grNo)
                                jj=jj+1;
                                these_ii_PLVs(jj)=ii;
                            end
                        end
                        
                        %Get delta phases
                        
                        mean_dp=[];
                        if ~isfield(handles_out,'drg')
                            Fs=20000;
                        else
                            Fs=handles_out.drg.session(1).draq_p.ActualRate;
                        end
                        
                        for ii_plv_per_mouse=1:length(these_ii_PLVs)
                            this_ii_PLV=these_ii_PLVs(ii_plv_per_mouse);
                            these_DPs=zeros(length(handles.drgbchoices.reference_electrodes)*length(handles.drgbchoices.other_electrodes),length(handles_out.PLV(this_ii_PLV).time));
                            ii_PLVs=0;
                            for ii_ref=handles.drgbchoices.reference_electrodes
                                for ii_oth=handles.drgbchoices.other_electrodes
                                    this_dp=zeros(1,length(handles_out.PLV(this_ii_PLV).time));
                                    if ii_ref<ii_oth
                                        this_dp(1,:)=handles_out.PLV(this_ii_PLV).delta_phase(:,ii_ref,ii_oth);
                                    else
                                        this_dp(1,:)=handles_out.PLV(this_ii_PLV).delta_phase(:,ii_oth,ii_ref);
                                    end
                                    ii_PLVs=ii_PLVs+1;
                                    these_DPs(ii_PLVs,:)=this_dp;
                                end
                            end
                            
                            these_DPs_dec=zeros(size(these_DPs,1),size(these_DPs,2)/no_time_points);
                            time_dec=zeros(1,size(these_DPs,2)/no_time_points);
                            time=handles_out.PLV(this_ii_PLV).time;
                            
                            for elec_pair=1:size(these_DPs,1)
                                for ii_dec_tp=1:size(these_DPs,2)/no_time_points
                                    these_DPs_dec(elec_pair,ii_dec_tp)=circ_mean(these_DPs(elec_pair,1+(ii_dec_tp-1)*no_time_points:ii_dec_tp*no_time_points)');
                                    if elec_pair==1
                                        time_dec(1,ii_dec_tp)=mean(time(1+(ii_dec_tp-1)*no_time_points:ii_dec_tp*no_time_points));
                                    end
                                end
                            end
                            
                            mean_dp(ii_plv_per_mouse,:)=circ_mean(these_DPs_dec);
                            
                        end
                        
                        if eventNo_ii==2
                            plot(time_dec,circ_mean(mean_dp), 'b');
                        else
                            plot(time_dec,circ_mean(mean_dp),'r');
                        end
                        
                        %                 ylim([0.2 1])
                        xlabel('Time (sec)')
                        ylabel('Delta phase')
                        title([handles.drgbchoices.bwlabels{bw_ii} ' ' handles.drgbchoices.percent_labels{per_ii}])
                        
                        %                 if (bw_ii==1)&(per_ii==size(handles.drgbchoices.percent_windows,1))
                        %                     text(0.95,-1.5,'S+','Color','r')
                        %                     text(0.9,-1.5,'S-','Color','b')
                        %                     text(0.9,-1.5,'p<0.05','Color','k')
                        %                 end
                        
                    end
                end
            end
            
            suptitle(['Delta phase for ' handles.drgbchoices.group_no_names{grNo}])
            
            
end






for grNo=1:3
    
    
    %Plot average PLV timecourse between brain regions
    
    
    
    
    
    for bw_ii=1:length(handles.drgbchoices.lowF)
        
        
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        
        
        hFig = figure(figNo);
        set(hFig, 'units','normalized','position',[.05 .05 .45 .7])
        
        this_sub=0;
        
        per_ii=1;
        
        for mouse_ii=1:max(handles.drgbchoices.mouse_no)
            
            PLV_exists=0;
            for ii=1:length(handles_out.PLV)
                if (handles_out.PLV(ii).per_ii==per_ii)&(handles_out.PLV(ii).bw_ii==bw_ii)...
                        &(handles_out.PLV(ii).group_no==grNo)&(handles_out.PLV(ii).mouse_no==mouse_ii)
                    PLV_exists=1;
                end
            end
            
            if PLV_exists==1
                this_sub=this_sub+1;
                subplot(4,2,this_sub)
                subplot(length(handles.drgbchoices.lowF),size(handles.drgbchoices.percent_windows,1),this_sub)
                hold on
                
                for eventNo_ii=length(handles.drgbchoices.evTypeNos):-1:1
                    
                    %Get the PLVs
                    this_ii_PLV=[];
                    
                    
                    for ii=1:length(handles_out.PLV)
                        if (handles_out.PLV(ii).eventNo_ii==eventNo_ii)&(handles_out.PLV(ii).per_ii==per_ii)&(handles_out.PLV(ii).bw_ii==bw_ii)...
                                &(handles_out.PLV(ii).group_no==grNo)&(handles_out.PLV(ii).mouse_no==mouse_ii)
                            this_ii_PLV=ii;
                        end
                    end
                    
                    if ~isempty(this_ii_PLV)
                        mean_plvs=[];
                        if ~isfield(handles_out,'drg')
                            Fs=20000;
                        else
                            Fs=handles_out.drg.session(1).draq_p.ActualRate;
                        end
                        
                        these_PLVs=zeros(length(handles.drgbchoices.reference_electrodes)*length(handles.drgbchoices.other_electrodes),length(handles_out.PLV(this_ii_PLV).time));
                        ii_PLVs=0;
                        for ii_ref=handles.drgbchoices.reference_electrodes
                            for ii_oth=handles.drgbchoices.other_electrodes
                                this_plv=zeros(1,length(handles_out.PLV(this_ii_PLV).time));
                                if ii_ref<ii_oth
                                    this_plv(1,:)=handles_out.PLV(this_ii_PLV).plv(:,ii_ref,ii_oth);
                                else
                                    this_plv(1,:)=handles_out.PLV(this_ii_PLV).plv(:,ii_oth,ii_ref);
                                end
                                ii_PLVs=ii_PLVs+1;
                                these_PLVs(ii_PLVs,:)=this_plv;
                            end
                        end
                        
                        %Calculate plv for 100 msec intervals
                        dt=0.1;
                        decimation_factor=40;
                        no_time_points=(Fs/decimation_factor)*dt;
                        these_PLVs_dec=zeros(size(these_PLVs,1),size(these_PLVs,2)/no_time_points);
                        time_dec=zeros(1,size(these_PLVs,2)/no_time_points);
                        time=handles_out.PLV(this_ii_PLV).time;
                        
                        for elec_pair=1:size(these_PLVs,1)
                            for ii_dec_tp=1:size(these_PLVs,2)/no_time_points
                                these_PLVs_dec(elec_pair,ii_dec_tp)=mean(these_PLVs(elec_pair,1+(ii_dec_tp-1)*no_time_points:ii_dec_tp*no_time_points));
                                if elec_pair==1
                                    time_dec(1,ii_dec_tp)=mean(time(1+(ii_dec_tp-1)*no_time_points:ii_dec_tp*no_time_points));
                                end
                            end
                        end
                        
                        mean_plv=mean(these_PLVs_dec);
                        
                        if eventNo_ii==2
                            plot(time_dec,mean_plv, 'b');
                        else
                            plot(time_dec,mean_plv, 'r');
                        end
                        
                        %Get PLS
                        these_PLVs=zeros(length(handles.drgbchoices.reference_electrodes),length(handles_out.PLV(this_ii_PLV).time));
                        ii_PLVs=0;
                        for ii_ref=1:handles.drgbchoices.no_electrodes
                            try
                                ii_oth=handles.drgbchoices.no_electrodes+1;
                                this_plv=zeros(1,length(handles_out.PLV(this_ii_PLV).time));
                                if ii_ref<ii_oth
                                    this_plv(1,:)=handles_out.PLV(this_ii_PLV).plv(:,ii_ref,ii_oth);
                                else
                                    this_plv(1,:)=handles_out.PLV(this_ii_PLV).plv(:,ii_oth,ii_ref);
                                end
                            catch
                            end
                            ii_PLVs=ii_PLVs+1;
                            these_PLVs(ii_PLVs,:)=this_plv;
                            
                        end
                        pls=prctile(these_PLVs(:),95);
                        plot([time_dec(1) time_dec(end)],[pls pls],'-k')
                        
                        ylim([0 1])
                        xlabel('Time (sec)')
                        ylabel('PLV')
                        title(['Mouse ' num2str(mouse_ii) ])
                        
                        %                         if (bw_ii==1)&(per_ii==size(handles.drgbchoices.percent_windows,1))
                        %                             text(0.95,-1.5,'S+','Color','r')
                        %                             text(0.9,-1.5,'S-','Color','b')
                        %                             text(0.9,-1.5,'p<0.05','Color','k')
                        %                         end
                    end
                end
            end
        end
        suptitle(['PLV for ' handles.drgbchoices.group_no_names{grNo} ' ' handles.drgbchoices.bwlabels{bw_ii}])
    end
end

pfft=1;
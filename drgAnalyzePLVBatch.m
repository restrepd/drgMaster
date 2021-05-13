%drgAnalyzePLVBatch analyzes the output from drgRunBatchPLVpar

close all
clear all

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

reference_window=[-1 0];
odor_window=[1.5 2.5];

[choiceFileName,choiceBatchPathName] = uigetfile({'drgbChoicesPLV*.m'},'Select the .m file with all the choices for analysis');

% [choiceFileName,choiceBatchPathName] = uigetfile();

[PLVFileName,PLVBatchPathName] = uigetfile({'*.mat'},'Select the .mat file with the results of drgRunBatchPLVpar');

load([PLVBatchPathName PLVFileName])

%Display the data


tempDirName=['temp' choiceFileName(12:end-2)];

addpath(choiceBatchPathName)
eval(['handles_in=' choiceFileName(1:end-2) ';'])

handles.drgbchoices=handles_in.drgbchoices;

figNo=0;

delta_PLV_odor_minus_ref=[];
delta_PLV_odor_minus_ref_per_mouse=[];
delta_phase_odor=[];
delta_phase_odor_per_mouse=[];
group_no_per_mouse=[];
which_mice=[];
no_mice_included=0;

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
                these_mouse_nos=[];
                jj=0;
                
                for ii=1:length(handles_out.PLV)
                    if (handles_out.PLV(ii).eventNo_ii==eventNo_ii)&(handles_out.PLV(ii).per_ii==per_ii)&(handles_out.PLV(ii).bw_ii==bw_ii)&(handles_out.PLV(ii).group_no==grNo)
                        jj=jj+1;
                        these_ii_PLVs(jj)=ii;
                        these_mouse_nos(jj)=handles_out.PLV(ii).mouse_no;
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
                
                CI = bootci(1000, {@mean, mean_plv})';
                CI(:,1)=mean(mean_plv)'-CI(:,1);
                CI(:,2)=CI(:,2)-mean(mean_plv)';
                
                
                if eventNo_ii==2
                    [hlCR, hpCR] = boundedline(time_dec,mean(mean_plv), CI, 'b');
                else
                    [hlCR, hpCR] = boundedline(time_dec,mean(mean_plv), CI, 'r');
                end
                
                %
                %                         if eventNo_ii==2
                %                             plot(time_dec,mean(mean_plv), 'b');
                %                         else
                %                             plot(time_dec,mean(mean_plv), 'r');
                %                         end
                
                %                         mean_mean_plv=mean(mean_plv);
                %                         delta_PLV_odor_minus_ref(grNo,bw_ii,per_ii,eventNo_ii)=mean(mean_mean_plv((time_dec>=odor_window(1))&(time_dec<=odor_window(2))))-...
                %                             mean(mean_mean_plv((time_dec>=reference_window(1))&(time_dec<=reference_window(2))));
                
                these_delta_PLV_odor_minus_ref_per_mouse=[];
                for mouse_ii=1:length(these_mouse_nos)
                    if no_mice_included==0
                        no_mice_included=1;
                        which_mice(1)=these_mouse_nos(mouse_ii);
                        this_mouse_ii=1;
                    else
                        this_mouse_ii=find(which_mice==these_mouse_nos(mouse_ii));
                        if isempty(this_mouse_ii)
                            no_mice_included=no_mice_included+1;
                            which_mice(no_mice_included)=these_mouse_nos(mouse_ii);
                            this_mouse_ii=no_mice_included;
                        end
                    end
                    this_mean_plv=zeros(size(mean_plv,2),1);
                    this_mean_plv(:,1)=mean_plv(mouse_ii,:);
                    %                             delta_PLV_odor_minus_ref_per_mouse(these_mouse_nos(mouse_ii),bw_ii,per_ii,eventNo_ii) = mean(this_mean_plv( (time_dec>=odor_window(1))&(time_dec<=odor_window(2)) ))-...
                    %                             mean( this_mean_plv((time_dec>=reference_window(1))&(time_dec<=reference_window(2))) );
                    delta_PLV_odor_minus_ref_per_mouse(this_mouse_ii,bw_ii,per_ii,eventNo_ii) = mean(this_mean_plv( (time_dec>=odor_window(1))&(time_dec<=odor_window(2)) )-...
                        this_mean_plv((time_dec>=reference_window(1))&(time_dec<=reference_window(2))) );
                    these_delta_PLV_odor_minus_ref_per_mouse(mouse_ii)=delta_PLV_odor_minus_ref_per_mouse(this_mouse_ii,bw_ii,per_ii,eventNo_ii);
                    group_no_per_mouse(this_mouse_ii)=grNo;
                end
                
                %                         delta_PLV_odor_minus_ref(grNo,bw_ii,per_ii,eventNo_ii)=mean(these_delta_PLV_odor_minus_ref_per_mouse);
                
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
                these_mouse_nos=[];
                jj=0;
                
                for ii=1:length(handles_out.PLV)
                    if (handles_out.PLV(ii).eventNo_ii==eventNo_ii)&(handles_out.PLV(ii).per_ii==per_ii)&(handles_out.PLV(ii).bw_ii==bw_ii)&(handles_out.PLV(ii).group_no==grNo)
                        jj=jj+1;
                        these_ii_PLVs(jj)=ii;
                        these_mouse_nos(jj)=handles_out.PLV(ii).mouse_no;
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
                
                CI = bootci(1000, {@circ_mean, mean_dp})';
                CI(:,1)=mean(mean_dp)'-CI(:,1);
                CI(:,2)=CI(:,2)-mean(mean_dp)';
                
                
                if eventNo_ii==2
                    [hlCR, hpCR] = boundedline(time_dec,circ_mean(mean_dp), CI, 'b');
                else
                    [hlCR, hpCR] = boundedline(time_dec,circ_mean(mean_dp), CI, 'r');
                end
                
                circ_mean_mean_dp=circ_mean(mean_dp);
                delta_phase_odor(grNo,bw_ii,per_ii,eventNo_ii)=circ_mean(circ_mean_mean_dp((time_dec>=odor_window(1))&(time_dec<=odor_window(2)))');
                
                
                for mouse_ii=1:length(these_mouse_nos)
                    
                    this_mouse_ii=find(which_mice==these_mouse_nos(mouse_ii));
                    
                    this_mean_dp=zeros(size(mean_dp,2),1);
                    this_mean_dp(:,1)=mean_dp(mouse_ii,:);
                    delta_phase_odor_per_mouse(this_mouse_ii,bw_ii,per_ii,eventNo_ii)=circ_mean(this_mean_dp((time_dec>=odor_window(1))&(time_dec<=odor_window(2)),1));
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
                        
                        %                         mean_plv=mean(these_PLVs_dec);
                        
                        %                         if eventNo_ii==2
                        %                             plot(time_dec,mean_plv, 'b');
                        %                         else
                        %                             plot(time_dec,mean_plv, 'r');
                        %                         end
                        
                        CI = bootci(1000, {@mean, these_PLVs_dec})';
                        CI(:,1)=mean(these_PLVs_dec)'-CI(:,1);
                        CI(:,2)=CI(:,2)-mean(these_PLVs_dec)';
                        
                        
                        if eventNo_ii==2
                            [hlCR, hpCR] = boundedline(time_dec,mean(these_PLVs_dec), CI, 'b');
                        else
                            [hlCR, hpCR] = boundedline(time_dec,mean(these_PLVs_dec), CI, 'r');
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

for grNo=1:3
    for bw_ii=1:length(handles.drgbchoices.lowF)
        for per_ii=size(handles.drgbchoices.percent_windows,1):-1:1
            for eventNo_ii=length(handles.drgbchoices.evTypeNos):-1:1
                delta_PLV_odor_minus_ref(grNo,bw_ii,per_ii,eventNo_ii)=mean(delta_PLV_odor_minus_ref_per_mouse(group_no_per_mouse==grNo,bw_ii,per_ii,eventNo_ii));
            end
        end
    end
end

drgbchoices=handles.drgbchoices;
save([PLVBatchPathName PLVFileName(1:end-4) '_out.mat'],'delta_PLV_odor_minus_ref','delta_phase_odor','drgbchoices','delta_PLV_odor_minus_ref_per_mouse'...
    ,'delta_phase_odor_per_mouse','group_no_per_mouse','which_mice');


%Now plot the  delta PLV for each odor pair/mouse
%(including all sessions for each mouse)
edges=[-0.6:0.05:0.6];
rand_offset=0.7;

for bw_ii=1:4    %for amplitude bandwidths (beta, low gamma, high gamma)
    
    glm_PLV=[];
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
    
    %             try
    %                 close(figNo+pacii)
    %             catch
    %             end
    %             hFig=figure(figNo+pacii);
    
    set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
    hold on
    
    bar_lab_loc=[];
    no_ev_labels=0;
    ii_gr_included=0;
    bar_offset = 0;
    
    for evNo=1:2
        
        for per_ii=2:-1:1
            
            for grNo=1:3
                bar_offset = bar_offset +1;
                
                %                         if sum(eventType==3)>0
                %                             bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(2-evNo);
                %                         else
                %                             bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(length(eventType)-evNo);
                %                         end
                %
                %                         these_offsets(per_ii)=bar_offset;
                bar_offset = bar_offset + 1;
                
                if (grNo==3)&(evNo==2)&(per_ii==1)
                    pffft=1;
                end
                
                %Get these PLV values
                %                 these_PLV=[];
                %                 ii_PLV=0;
                %                 for ii=1:length(FileName)
                %                     these_PLVs=zeros(1,sum(group_no_per_mouse==grNo));
                %                     these_PLVs(1,:)=all_files(ii).delta_PLV_odor_minus_ref_per_mouse(group_no_per_mouse==grNo,bwii,per_ii,evNo);
                %                     these_PLV(ii_PLV+1:ii_PLV+length(these_PLVs))=these_PLVs;
                %                     ii_PLV=ii_PLV+length(these_PLVs);
                %                 end
                
                
                if (grNo==3)&(evNo==2)&(per_ii==1)
                    pffft=1;
                end
                
                
                this_delta_PLV=delta_PLV_odor_minus_ref(grNo,bw_ii,per_ii,evNo);
                
                switch grNo
                    case 1
                        bar(bar_offset,mean(this_delta_PLV),'g','LineWidth', 3,'EdgeColor','none')
                    case 2
                        bar(bar_offset,mean(this_delta_PLV),'b','LineWidth', 3,'EdgeColor','none')
                    case 3
                        bar(bar_offset,mean(this_delta_PLV),'y','LineWidth', 3,'EdgeColor','none')
                end
                
                these_delta_PLVs=delta_PLV_odor_minus_ref_per_mouse(group_no_per_mouse==grNo,bw_ii,per_ii,evNo);
                plot(bar_offset,mean(these_delta_PLVs),'ok','MarkerSize',10)
                
                %Violin plot
                
                %                 [mean_out, CIout]=drgViolinPoint(these_PLV,edges,bar_offset,rand_offset,'k','k',3);
                %                 CI = bootci(1000, {@mean, these_PLV},'type','cper');
                %                 plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                plot(bar_offset*ones(1,length(these_delta_PLVs)),these_delta_PLVs,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
                %
                
                %                                 %Save data for glm and ranksum
                
                glm_PLV.data(glm_ii+1:glm_ii+length(these_delta_PLVs))=these_delta_PLVs;
                glm_PLV.group(glm_ii+1:glm_ii+length(these_delta_PLVs))=grNo*ones(1,length(these_delta_PLVs));
                glm_PLV.perCorr(glm_ii+1:glm_ii+length(these_delta_PLVs))=per_ii*ones(1,length(these_delta_PLVs));
                glm_PLV.event(glm_ii+1:glm_ii+length(these_delta_PLVs))=evNo*ones(1,length(these_delta_PLVs));
                glm_ii=glm_ii+length(these_delta_PLVs);
                
                id_ii=id_ii+1;
                input_data(id_ii).data=these_delta_PLVs;
                input_data(id_ii).description=[group_legend{grNo} ' ' evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
                
                
            end
            bar_offset = bar_offset + 2;
            
        end
        bar_offset = bar_offset + 3;
        
    end
    
    title(['Average delta PLV for each odor_pair/mouse for ' bandwidth_names{bw_ii}])
    
    
    %Annotations identifying groups
    x_interval=0.8/ii_gr_included;
    for ii=1:ii_gr_included
        annotation('textbox',[0.7*x_interval+x_interval*(ii-1) 0.7 0.3 0.1],'String',handles_drgb.drgbchoices.group_no_names{ groups_included(ii)},'FitBoxToText','on');
    end
    
    %Proficient/Naive annotations
    annotation('textbox',[0.15 0.8 0.3 0.1],'String','Proficient','FitBoxToText','on','Color','r','LineStyle','none');
    annotation('textbox',[0.15 0.75 0.3 0.1],'String','Naive','FitBoxToText','on','Color','b','LineStyle','none');
    
    
    xticks([2 4 6 10 12 14 21 23 25 29 31 33])
    xticklabels({'nwtS+', 'nHS+', 'nKOS+', 'pwtS+', 'pHS+', 'pKOS+', 'nwtS-', 'nHS-', 'nKOS-', 'pwtS-', 'pHS-', 'pKOS-'})
    
    ylabel('delta PLV')
    
    
    %Perform the glm
    fprintf(1, ['\n\nglm for delta PLV for each odor pair/mouse  for' bandwidth_names{bw_ii} '\n'])
    tbl = table(glm_PLV.data',glm_PLV.group',glm_PLV.perCorr',glm_PLV.event',...
        'VariableNames',{'MI','group','perCorr','event'});
    mdl = fitglm(tbl,'MI~group+perCorr+event+perCorr*group*event'...
        ,'CategoricalVars',[2,3,4])
    
    
    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for delta PLV for each odor pair for ' bandwidth_names{bw_ii} ' \n'])
    [output_data] = drgMutiRanksumorTtest(input_data);
    
    
end

fprintf(1, ['\ndrgAnalyzePLVBatcg run for ' PLVFileName '\n\n']);
 
pfft=1;
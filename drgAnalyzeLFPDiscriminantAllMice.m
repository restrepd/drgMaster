function drgAnalyzeLFPDiscriminantAllMice
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

handles_out=[];

[fname,pname,nCancel] = uigetfile({'Discriminant_*.mat'},'Select the per mouse discriminant output file ...');
if nCancel
    inputPath = [pname,fname];
    pnameStart = pname;
    %save('timeCorr_cfg.mat','pnameStart','-append');
else
    error('Cancelled')
end

discriminant_name=[pname fname];
load(discriminant_name)

%You can choose to process only a subset of mice
mice_included=[1:length(handles_out.discriminant_PACwavepower)];

%Used to troubleshoot the IAMO JL
% mice_included=[1 3 4];

figNo=0;
%Plot average percent correct for the LDA for peak and trough for
%wavelet power referenced to PAC phase
t=handles_out.t_power;

for PACii=1:length(handles_out.drgbchoices.PACburstLowF)
    
    for percent_correct_ii=1:2
        
        for groupNo=1:max(handles_out.drgbchoices.group_no)
            

            p_val_lick=handles_all.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).p_val_lick;
            p_val_trough=handles_all.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).p_val_trough;
            p_val_peak=handles_all.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).p_val_peak;
            

            
            %Determine decision times for the all mouse data
            
            
            wing=1;
            
            jj_start=find(t>=t_odor_arrival,1,'first');
            %Find the discrimination time
            
            
            found_disc_t=0;
            while (found_disc_t==0)&(jj_start<length(t))
                ii_next=find(p_val_lick(1,jj_start:end)<=0.05,1,'first');
                if isempty(ii_next)
                    jj_start=length(t);
                    ii_next=1;
                else
                    if jj_start+ii_next+wing>length(t)
                        found_disc_t=1;
                    else
                        if sum(p_val_lick(1,jj_start+ii_next:jj_start+ii_next+wing)>0.05)>=1
                            jj_start=jj_start+ii_next;
                        else
                            found_disc_t=1;
                        end
                    end
                end
            end
            t_detect_licks_all=t(jj_start+ii_next-1)-t_odor_arrival;
            
            
            fprintf(1, ['Lick discrimination time for all mouse analysis (sec) %d\n'],mean(t_detect_licks_all))
            
            
            
            jj_start=find(t>=t_odor_arrival,1,'first');
            %Find the discrimination time
            
            found_disc_t=0;
            while (found_disc_t==0)&(jj_start<length(t))
                ii_next=find(p_val_trough(1,jj_start:end)<=0.05,1,'first');
                if isempty(ii_next)
                    jj_start=length(t);
                    ii_next=1;
                else
                    if jj_start+ii_next+wing>length(t)
                        found_disc_t=1;
                    else
                        if sum(p_val_trough(1,jj_start+ii_next:jj_start+ii_next+wing)>0.05)>=1
                            jj_start=jj_start+ii_next;
                        else
                            found_disc_t=1;
                        end
                    end
                end
            end
            t_detect_trough_all=t(jj_start+ii_next-1)-t_odor_arrival;
            
            fprintf(1, ['Trough discrimination time for all mouse analysis (sec) %d\n'],mean(t_detect_trough_all))
            
            
            jj_start=find(t>=t_odor_arrival,1,'first');
            
            %Find the discrimination time
            
            
            found_disc_t=0;
            while (found_disc_t==0)&(jj_start<length(t))
                ii_next=find(p_val_peak(1,jj_start:end)<=0.05,1,'first');
                if isempty(ii_next)
                    jj_start=length(t);
                    ii_next=1;
                else
                    if jj_start+ii_next+wing>length(t)
                        found_disc_t=1;
                    else
                        if sum(p_val_peak(1,jj_start+ii_next:jj_start+ii_next+wing)>0.05)>=1
                            jj_start=jj_start+ii_next;
                        else
                            found_disc_t=1;
                        end
                    end
                end
            end
            t_detect_peak_all=t(jj_start+ii_next-1)-t_odor_arrival;
           
            fprintf(1, ['Peak discrimination time for all mouse analysis (sec) %d\n'],t_detect_peak_all)
            
           %Now do the per mouse analysis
            
            %Gather all the data
            no_mice=0;
            no_mice_included=0;
            all_discriminant_p_val_lick=zeros(max(handles_out.drgbchoices.mouse_no),length(t));
            all_discriminant_p_val_peak=zeros(max(handles_out.drgbchoices.mouse_no),length(t));
            all_discriminant_p_val_trough=zeros(max(handles_out.drgbchoices.mouse_no),length(t));
            
            for mouseNo=1:length(handles_out.discriminant_PACwavepower)
                if sum(mice_included==mouseNo)>0
                    try
                        if handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated==1
                            per_ii=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).no_trials;
                            no_mice=no_mice+1;
                            if (per_ii>=20)&(sum(no_mice==mice_excluded)==0)
                                no_mice_included=no_mice_included+1;
                                all_discriminant_p_val_lick(no_mice_included,:)=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).p_val_lick;
                                all_discriminant_p_val_peak(no_mice_included,:)=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).p_val_peak;
                                all_discriminant_p_val_trough(no_mice_included,:)=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).p_val_trough;
                            end
                        end
                    catch
                    end
                end
            end
            
         
            no_mice_per(percent_correct_ii)=no_mice_included;
            

            wing=1;
            t_detect=zeros(1,no_mice_included);
            jj_start=find(t>=t_odor_arrival,1,'first');
            %Find the discrimination time
            for ii=1:no_mice_included
                this_p_val_licks=[];
                this_p_val_licks=all_discriminant_p_val_lick(ii,:);
                found_disc_t=0;
                while (found_disc_t==0)&(jj_start<length(t))
                    ii_next=find(this_p_val_licks(1,jj_start:end)<=0.05,1,'first');
                    if isempty(ii_next)
                        jj_start=length(t);
                        ii_next=1;
                    else
                        if jj_start+ii_next+wing>length(t)
                            found_disc_t=1;
                        else
                            if sum(this_p_val_licks(1,jj_start+ii_next:jj_start+ii_next+wing)>0.05)>=1
                                jj_start=jj_start+ii_next;
                            else
                                found_disc_t=1;
                            end
                        end
                    end
                end
                t_detect(ii)=t(jj_start+ii_next-1)-t_odor_arrival;
            end
            
            fprintf(1, ['Mean lick discrimination time (sec) %d\n'],mean(t_detect))
            
            lick_t_detect=t_detect;
            
            
            
            t_detect=zeros(1,no_mice_included);
            jj_start=find(t>=t_odor_arrival,1,'first');
            %Find the discrimination time
            for ii=1:no_mice_included
                this_p_val_troughs=[];
                this_p_val_troughs=all_discriminant_p_val_trough(ii,:);
                found_disc_t=0;
                while (found_disc_t==0)&(jj_start<length(t))
                    ii_next=find(this_p_val_troughs(1,jj_start:end)<=0.05,1,'first');
                    if isempty(ii_next)
                        jj_start=length(t);
                        ii_next=1;
                    else
                        if jj_start+ii_next+wing>length(t)
                            found_disc_t=1;
                        else
                            if sum(this_p_val_troughs(1,jj_start+ii_next:jj_start+ii_next+wing)>0.05)>=1
                                jj_start=jj_start+ii_next;
                            else
                                found_disc_t=1;
                            end
                        end
                    end
                end
                t_detect(ii)=t(jj_start+ii_next-1)-t_odor_arrival;
            end
            
            trough_t_detect=t_detect;
            
            fprintf(1, ['Trough discrimination time (sec) %d\n'],mean(t_detect))
            
            
            
            t_detect=zeros(1,no_mice_included);
            jj_start=find(t>=t_odor_arrival,1,'first');
            %Find the discrimination time
            for ii=1:no_mice_included
                this_p_val_peaks=[];
                this_p_val_peaks=all_discriminant_p_val_peak(ii,:);
                found_disc_t=0;
                while (found_disc_t==0)&(jj_start<length(t))
                    ii_next=find(this_p_val_peaks(1,jj_start:end)<=0.05,1,'first');
                    if isempty(ii_next)
                        jj_start=length(t);
                        ii_next=1;
                    else
                        if jj_start+ii_next+wing>length(t)
                            found_disc_t=1;
                        else
                            if sum(this_p_val_peaks(1,jj_start+ii_next:jj_start+ii_next+wing)>0.05)>=1
                                jj_start=jj_start+ii_next;
                            else
                                found_disc_t=1;
                            end
                        end
                    end
                end
                t_detect(ii)=t(jj_start+ii_next-1)-t_odor_arrival;
            end
            
            fprintf(1, ['Peak discrimination time (sec) %d\n'],mean(t_detect))
            
            peak_t_detect=t_detect;
            
            
            
            %Plot decision time reationship
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            hFig=figure(figNo);
            hold on
            
            %Licks
            mean_lick_t_detect=mean(lick_t_detect)';
            handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).mean_lick_t_detect=mean_lick_t_detect;
            CIlicks = bootci(1000, {@mean, lick_t_detect})';
            p1=bar(1,mean_lick_t_detect,'EdgeColor','k','FaceColor',[0.7 0.7 0.7]);
            plot(1*ones(1,length(lick_t_detect)),lick_t_detect,'ok')
            plot([1 1],CIlicks,'-k')
            plot(1,t_detect_licks_all,'ok','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
            handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).t_detect_licks_all=t_detect_licks_all;
            
            %Trough
            mean_trough_t_detect=mean(trough_t_detect)';
            handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).mean_trough_t_detect=mean_trough_t_detect;
            CItrough = bootci(1000, {@mean, trough_t_detect})';
            p2=bar(2,mean_trough_t_detect,'b');
            plot(2*ones(1,length(trough_t_detect)),trough_t_detect,'ok')
            plot([2 2],CItrough,'-k')
            plot(2,t_detect_trough_all,'ok','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
            handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).t_detect_trough_all=t_detect_trough_all;
            
            %Peak
            mean_peak_t_detect=mean(peak_t_detect)';
            handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).mean_peak_t_detect=mean_peak_t_detect;
            CIpeak = bootci(1000, {@mean, peak_t_detect})';
            p3=bar(3,mean_peak_t_detect,'r');
            plot(3*ones(1,length(peak_t_detect)),peak_t_detect,'ok')
            plot([3 3],CIpeak,'-k')
            plot(3,t_detect_peak_all,'ok','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
            handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).t_detect_peak_all=t_detect_peak_all;
            
            title(['Decision times for Theta/' handles_out.drgbchoices.PACnames{PACii} ' ' handles_out.drgbchoices.group_no_names{groupNo}  ' ' handles_out.drgbchoices.per_lab{percent_correct_ii}])
            
            ylabel('Decision time (sec)')
            legend([p1 p2 p3],{'Lick','Peak','Trough'})
            
            if (PACii==3)&(percent_correct_ii==1)&(groupNo==1)
                pffft=1;
            end
            pffft=1;
        end
    end
end


%Now perform the analysis for the area under the curve
auc_from=0.1;
auc_to=2.5;
no_wins=2;
window_start=[-1 0.5];
window_end=[0 2.5];

for PACii=1:length(handles_out.drgbchoices.PACburstLowF)
    p_correct_stats=[];
    ii_stats=0;
    p_dim_stats=[];
    ii_dim_stats=0;
    glm_ii=0;
    glm_correct=[];
    glm_dim_ii=0;
    glm_dim=[];
    for percent_correct_ii=1:2
        
        for groupNo=1:max(handles_out.drgbchoices.group_no)
            
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
                            all_dimensionality_peak(no_mice_included,:)=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).dimensionality_peak;
                            all_dimensionality_trough(no_mice_included,:)=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).dimensionality_trough;
                        end
                    end
                catch
                end
            end
            
            
            no_mice_per(percent_correct_ii)=no_mice_included;
            
            auc_peak_shuffled=[];
            for mouseNo=1:no_mice_included
                auc_peak_shuffled=[auc_peak_shuffled (mean(all_discriminant_correct_shuffled_peak(mouseNo,(t>=auc_from)&(t<=auc_to)),2)-50)/50];
            end
            
            auc_trough_shuffled=[];
            for mouseNo=1:no_mice_included
                auc_trough_shuffled=[auc_trough_shuffled (mean(all_discriminant_correct_shuffled_trough(mouseNo,(t>=auc_from)&(t<=auc_to)),2)-50)/50];
            end
            
            auc_trough=[];
            for mouseNo=1:no_mice_included
                auc_trough=[auc_trough (mean(all_discriminant_correct_trough(mouseNo,(t>=auc_from)&(t<=auc_to)),2)-50)/50];
            end
            
            auc_peak=[];
            for mouseNo=1:no_mice_included
                auc_peak=[auc_peak (mean(all_discriminant_correct_peak(mouseNo,(t>=auc_from)&(t<=auc_to)),2)-50)/50];
            end
           
            
            %Now get the all mouse auc values
            
            all_mouse_discriminant_correct_peak=handles_all.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_peak;
            auc_peak_all=(mean(all_mouse_discriminant_correct_peak(1,(t>=auc_from)&(t<=auc_to)),2)-50)/50;
            all_mouse_discriminant_correct_trough=handles_all.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_trough;
            auc_trough_all=(mean(all_mouse_discriminant_correct_trough(1,(t>=auc_from)&(t<=auc_to)),2)-50)/50;
            all_mouse_discriminant_correct_shuffled_peak=handles_all.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_shuffled_peak;
            auc_peak_shuffled_all=(mean(all_mouse_discriminant_correct_shuffled_peak(1,(t>=auc_from)&(t<=auc_to)),2)-50)/50;
            all_mouse_discriminant_correct_shuffled_trough=handles_all.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_shuffled_trough;
            auc_trough_shuffled_all=(mean(all_mouse_discriminant_correct_shuffled_trough(1,(t>=auc_from)&(t<=auc_to)),2)-50)/50;
            
            %Plot the bar graph for auc
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            hFig=figure(figNo);
            hold on
            
            %Trough
            mean_trough_auc=mean(auc_trough)';
            handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).mean_trough_auc=mean_trough_auc;
            CItrough = bootci(1000, {@mean, auc_trough})';
            p1=bar(2,mean_trough_auc,'b');
            plot(2*ones(1,length(auc_trough)),auc_trough,'ok')
            plot([2 2],CItrough,'-k')
            plot(2,auc_trough_all,'ok','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
            handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).auc_trough_all=auc_trough_all;
            
            %Peak
            mean_peak_auc=mean(auc_peak)';
            handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).mean_peak_auc=mean_peak_auc;
            CIpeak = bootci(1000, {@mean, auc_peak})';
            p2=bar(3,mean_peak_auc,'r');
            plot(3*ones(1,length(auc_peak)),auc_peak,'ok')
            plot([3 3],CIpeak,'-k')
            plot(3,auc_peak_all,'ok','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
            handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).auc_peak_all=auc_peak_all;
            
            %Trough shuffled
            mean_trough_shuffled_auc=mean(auc_trough_shuffled)';
            handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).mean_trough_shuffled_auc=mean_trough_shuffled_auc;
            CItrough_shuffled = bootci(1000, {@mean, auc_trough_shuffled})';
            p3=bar(5,mean_trough_shuffled_auc,'EdgeColor',[0.7 0.7 1],'FaceColor',[0.7 0.7 1]);
            plot(5*ones(1,length(auc_trough_shuffled)),auc_trough_shuffled,'ok')
            plot([5 5],CItrough_shuffled,'-k')
            plot(5,auc_trough_shuffled_all,'ok','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
            handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).auc_trough_shuffled_all=auc_trough_shuffled_all;
            
            %Peak
            mean_peak_shuffled_auc=mean(auc_peak_shuffled)';
            handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).mean_peak_shuffled_auc=mean_peak_shuffled_auc;
            CIpeak_shuffled = bootci(1000, {@mean, auc_peak_shuffled})';
            p4=bar(6,mean_peak_shuffled_auc,'EdgeColor',[1 0.7 0.7],'FaceColor',[1 0.7 0.7]);
            plot(6*ones(1,length(auc_peak_shuffled)),auc_peak_shuffled,'ok')
            plot([6 6],CIpeak_shuffled,'-k')
            plot(6,auc_peak_shuffled_all,'ok','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
            handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).auc_peak_shuffled_all=auc_peak_shuffled_all;
            
            title(['Area under the curve for Theta/' handles_out.drgbchoices.PACnames{PACii} ' ' handles_out.drgbchoices.group_no_names{groupNo}  ' ' handles_out.drgbchoices.per_lab{percent_correct_ii}])
            
            ylabel('AUC')
            legend([p1 p2 p3 p4],{'Trough','Peak','Trough shuffled','Peak shuffled'})
            ylim([-0.1 1.1])
            
            if (PACii==3)&(percent_correct_ii==1)&(groupNo==1)
                pffft=1;
            end
        end
    end
    
end

save([pname 'LDA_all_m' fname(13:end)],'handles_outp','-v7.3')

 fprintf(1, ['Processed drgAnalyzeLFPDiscriminantAllMice with the following input files\n'])
 
 fprintf(1, discriminant_name_all)
 
 fprintf(1, discriminant_name_all)
           








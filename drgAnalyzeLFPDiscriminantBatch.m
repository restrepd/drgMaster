function drgAnalyzeLFPDiscriminantBatch

[fname,pname,nCancel] = uigetfile({'Discriminant_*.mat'},'Select the perceptron LFP batch output file ...');
if nCancel
    inputPath = [pname,fname];
    pnameStart = pname;
    %save('timeCorr_cfg.mat','pnameStart','-append');
else
    error('Cancelled')
end

discriminant_name=[pname fname];
load(discriminant_name)

figNo=0;

t=handles_out.t;

%Gather the data to be plotted


for groupNo=1:max(handles_out.drgbchoices.group_no)
    figNo=figNo+1
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    
    set(hFig, 'units','normalized','position',[.2 .2 .7 .7])
    
     
    ii_plot=0;
    for percent_correct_ii=1:2
        for bwii=1:4
            ii_plot=ii_plot+1;
            subplot(2,4,ii_plot);
            hold on
            no_mice=0;
            all_discriminant_correct=zeros(max(handles_out.drgbchoices.mouse_no),length(t));
            all_discriminant_correct_shuffled=zeros(max(handles_out.drgbchoices.mouse_no),length(t));
            for mouseNo=1:max(handles_out.drgbchoices.mouse_no)
                if isfield(handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii),'discriminant_calulated')
                    handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated=...
                        handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calulated;
                end
                if handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated==1
                    per_ii=handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).no_trials;
                    if per_ii>=20
                        no_mice=no_mice+1;
                        all_discriminant_correct(no_mice,:)=handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).discriminant_correct;
                        all_discriminant_correct_shuffled(no_mice,:)=handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).discriminant_correct_shuffled;
                    end
                end
            end
             
            no_mice_per(percent_correct_ii)=no_mice;
            
            all_discriminant_correct_shuffled=all_discriminant_correct_shuffled(1:no_mice,:);
            mean_dcsh=mean(all_discriminant_correct_shuffled,1)';
            CIdcsh = bootci(1000, {@mean, all_discriminant_correct_shuffled})';
            CIdcsh(:,1)=mean_dcsh-CIdcsh(:,1);
            CIdcsh(:,2)=CIdcsh(:,2)-mean_dcsh;
            [hlCR, hpCR] = boundedline(t,mean_dcsh, CIdcsh, 'b');
            
            all_discriminant_correct=all_discriminant_correct(1:no_mice,:);
            mean_dc=mean(all_discriminant_correct,1)';
            CIdc = bootci(1000, {@mean, all_discriminant_correct})';
            CIdc(:,1)=mean_dc-CIdc(:,1);
            CIdc(:,2)=CIdc(:,2)-mean_dc;
            [hlCR, hpCR] = boundedline(t,mean_dc, CIdc, 'r');
            
            
            %Odor on markers
            plot([0 0],[0 100],'-k')
            odorhl=plot([0 2.5],[20 20],'-k','LineWidth',5);
            plot([2.5 2.5],[0 100],'-k')
            
            title([handles_out.drgbchoices.bwlabels{bwii} ])
            
            xlabel('Time (sec)')
            ylabel(['% correct '  handles_out.drgbchoices.per_lab{percent_correct_ii}])
        end
    end
    suptitle(['All trials. Group: ' handles_out.drgbchoices.group_no_names{groupNo} ' # of mice: ' num2str(no_mice_per(1)) ' ' handles_out.drgbchoices.per_lab{1} ' ' num2str(no_mice_per(2)) ' ' handles_out.drgbchoices.per_lab{2}])
end



%Now do the analysis for each event

%Which tested event do all events belong to? For example Hit and Miss
%belong to S+
compare_to_event=[5 5 5 11 11 11];

if ~isfield(handles_out,'drg')
    handles_out.drg.draq_d.eventlabels{3}='Hit';
    handles_out.drg.draq_d.eventlabels{5}='S+';
    handles_out.drg.draq_d.eventlabels{7}='Miss';
    handles_out.drg.draq_d.eventlabels{9}='CR';
    handles_out.drg.draq_d.eventlabels{11}='S-';
    handles_out.drg.draq_d.eventlabels{13}='FA';
end
    
    

for evNo=1:length(handles_out.drgbchoices.evTypeNos)
    compare_event_ii=find(handles_out.drgbchoices.events_to_discriminate==compare_to_event(evNo));
    
        for groupNo=1:max(handles_out.drgbchoices.group_no)
            figNo=figNo+1
            try
                close(figNo)
            catch
            end
            hFig=figure(figNo);
            
            
            set(hFig, 'units','normalized','position',[.2 .2 .7 .7])
            
             
            ii_plot=0;
            for percent_correct_ii=1:2
                for bwii=1:4
                    ii_plot=ii_plot+1;
                    subplot(2,4,ii_plot);
                    hold on
                    
                    decoding_targets=zeros(10000,1);
                    test_out=zeros(10000,length(t));
                    no_trials=0;
                    
                    for mouseNo=1:max(handles_out.drgbchoices.mouse_no)
                        
                        if handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated==1
                            per_ii=handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).no_trials;
                            if per_ii>=20
                                these_events=handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).which_events(evNo,:);
                                if bwii==1
                                  fprintf(1, ['For ' handles_out.drg.draq_d.eventlabels{handles_out.drgbchoices.evTypeNos(evNo)} ', mouse no %d ' handles_out.drgbchoices.group_no_names{groupNo} ' ' handles_out.drgbchoices.per_lab{percent_correct_ii} ' there are %d trials\n'], mouseNo, sum(these_events));  
                                end
                                decoding_targets(no_trials+1:no_trials+sum(these_events),1)=handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).per_targets(compare_event_ii,logical(these_events));
                                for tt=1:length(t)
                                    test_out(no_trials+1:no_trials+sum(these_events),tt)=handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).test_out_per_timepoint(compare_event_ii,logical(these_events),tt);
                                end
                                no_trials=no_trials+sum(these_events);
                            end
                        end
                    end
                    
                    percent_correct=zeros(1,length(t));
                    for tt=1:length(t)
                        percent_correct(tt)=100*sum(test_out(1:no_trials,tt).*decoding_targets(1:no_trials,1))/sum(decoding_targets);
                    end
                    plot(t,percent_correct,'-r')
                    
                    %Try to bootstrap a CI here
                    
                    
                    %Odor on markers
                    plot([0 0],[0 100],'-k')
                    odorhl=plot([0 2.5],[20 20],'-k','LineWidth',5);
                    plot([2.5 2.5],[0 100],'-k')
                    
                    title([handles_out.drgbchoices.bwlabels{bwii} ])
                    
                    xlabel('Time (sec)')
                    ylabel(['% correct '  handles_out.drgbchoices.per_lab{percent_correct_ii}])
                end
            end
          
            suptitle([handles_out.drg.draq_d.eventlabels{handles_out.drgbchoices.evTypeNos(evNo)} ' Group: ' handles_out.drgbchoices.group_no_names{groupNo} ' # of mice: ' num2str(no_mice_per(1)) ' ' handles_out.drgbchoices.per_lab{1} ' ' num2str(no_mice_per(2)) ' ' handles_out.drgbchoices.per_lab{2}])
          
        end
    
    
end

%Calculate lick rate and compare to percent correct
pfft=1


for evNo=1:length(handles_out.drgbchoices.evTypeNos)
    compare_event_ii=find(handles_out.drgbchoices.events_to_discriminate==compare_to_event(evNo));
    
        for groupNo=1:max(handles_out.drgbchoices.group_no)
            figNo=figNo+1
            try
                close(figNo)
            catch
            end
            hFig=figure(figNo);
            
            
            set(hFig, 'units','normalized','position',[.2 .2 .7 .7])
            
             
            ii_plot=0;
            for percent_correct_ii=1:2
                for bwii=1:4
                    ii_plot=ii_plot+1;
                    subplot(2,4,ii_plot);
                    hold on
                    
                    lick_rate_per_trial=zeros(10000,length(t));
                    no_trials=0;
                    
                    for mouseNo=1:max(handles_out.drgbchoices.mouse_no)
                        
                        if handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated==1
                            per_ii=handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).no_trials;
                            if per_ii>=20
                                these_events=handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).which_events(evNo,:);
                                for events=1:length(these_events)
                                   if these_events(events)==1
                                       no_trials=no_trials+1;
                                       for lick_ii=1:handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).stamped_lick_ii(events)
                                           this_lick_t= handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).stamped_lick_times(events,lick_ii);
                                           dt=t(2)-t(1);
                                           if (this_lick_t>=t(1)-dt/2)&(this_lick_t<=t(end)+dt/2)
                                               tt=ceil((this_lick_t-t(1))/dt);
                                              if (tt>=1)&(tt<=length(t))
                                               lick_rate_per_trial(no_trials,tt)=lick_rate_per_trial(no_trials,tt)+(1/dt);
                                              end
                                           end
                                       end
                                   end
                                end
                              
                            end
                        end
                    end
                    
                    try
                        mean_lick_rate=mean(lick_rate_per_trial(1:no_trials,:),1)';
                        CImlr = bootci(1000, {@mean, lick_rate_per_trial(1:no_trials,:)})';
                        CImlr(:,1)=mean_lick_rate-CImlr(:,1);
                        CImlr(:,2)=CImlr(:,2)-mean_lick_rate;
                        [hlCR, hpCR] = boundedline(t,mean_lick_rate, CImlr, 'r');
                        
                    catch
                    end
                    %Odor on markers
                    plot([0 0],[0 10],'-k')
                    odorhl=plot([0 2.5],[1 1],'-k','LineWidth',5);
                    plot([2.5 2.5],[0 10],'-k')
                    
                    title([handles_out.drgbchoices.bwlabels{bwii} ])
                    
                    xlabel('Time (sec)')
                    ylabel(['Lick rate (Hz) '  handles_out.drgbchoices.per_lab{percent_correct_ii}])
                    
                end
            end
          
            suptitle([handles_out.drg.draq_d.eventlabels{handles_out.drgbchoices.evTypeNos(evNo)} ' Group: ' handles_out.drgbchoices.group_no_names{groupNo} ' # of mice: ' num2str(no_mice_per(1)) ' ' handles_out.drgbchoices.per_lab{1} ' ' num2str(no_mice_per(2)) ' ' handles_out.drgbchoices.per_lab{2}])
          
        end
    
    
end

pffft=1;


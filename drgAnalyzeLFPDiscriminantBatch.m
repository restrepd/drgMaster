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
                 
                if handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calulated==1
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
    suptitle(['Group: ' handles_out.drgbchoices.group_no_names{groupNo} ' # of mice: ' num2str(no_mice_per(1)) ' ' handles_out.drgbchoices.per_lab{1} ' ' num2str(no_mice_per(2)) ' ' handles_out.drgbchoices.per_lab{2}])
end
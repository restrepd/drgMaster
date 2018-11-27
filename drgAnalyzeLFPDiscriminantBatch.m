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

for mouseNo=1:max(handles_out.drgbchoices.mouse_no)
    for groupNo=1:max(handles_out.drgbchoices.group_no)
        for percent_correct_ii=1:2
            for bwii=1:length(handles_out.drgbchoices.lowF)
                if handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calulated==1
                    
                        per_ii=handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).no_trials;
                        
                        if per_ii>=20
                            figNo=figNo+1
                            try
                                close(figNo)
                            catch
                            end
                            
                            figure(figNo)
                            
                            hold on
                            
                            discriminant_correct_shuffled=zeros(1,length(t));
                            discriminant_correct_shuffled(1,:)=handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).discriminant_correct_shuffled;
                            discriminant_correct=zeros(1,length(t));
                            discriminant_correct(1,:)=handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).discriminant_correct;
                            
                            
                            per95=prctile(discriminant_correct_shuffled,95);
                            per5=prctile(discriminant_correct_shuffled,5);
                            CIsh=[mean(discriminant_correct_shuffled)-per5 per95-mean(discriminant_correct_shuffled)]';
                            [hlCR, hpCR] = boundedline([t(1) t(end)],[mean(discriminant_correct_shuffled) mean(discriminant_correct_shuffled)], CIsh', 'r');
                            
                            plot(t',discriminant_correct,'-k')
                            
                            %Odor on markers
                            plot([0 0],[0 100],'-k')
                            odorhl=plot([0 2.5],[20 20],'-k','LineWidth',5);
                            plot([2.5 2.5],[0 100],'-k')
                            
                            title(['% correct for ' handles_out.drgbchoices.bwlabels{bwii} ' mouse No ' num2str(mouseNo) ' ' handles_out.drgbchoices.per_lab{percent_correct_ii} ' ' handles_out.drgbchoices.group_no_names{groupNo}])
                       
                            xlabel('Time (sec)')
                            ylabel('Percent correct')
                        else
                            fprintf(1, ['Discriminant not processed for mouse No %d ' handles_out.drgbchoices.per_lab{percent_correct_ii} ' because there were only %d trials (fewer than 20 trials)\n'],mouseNo,per_ii);
                        end
                   
                end
            end
        end
    end
end
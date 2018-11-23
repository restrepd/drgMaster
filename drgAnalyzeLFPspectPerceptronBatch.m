function drgAnalyzeLFPspectPerceptronBatch

[fname,pname,nCancel] = uigetfile({'LFP_per_*.mat'},'Select the perceptron LFP batch output file ...');
if nCancel
    inputPath = [pname,fname];
    pnameStart = pname;
    %save('timeCorr_cfg.mat','pnameStart','-append');
else
    error('Cancelled')
end

per_name=[pname fname];
load(per_name)

figNo=0;

t=handles_out.t;

for mouseNo=1:handles_out.no_mice_with_data
    
    for percent_correct_ii=1:2
        
        if ~isempty(handles_out.perceptron_per_mouse(mouseNo).percent_correct)
            if percent_correct_ii<=length(handles_out.perceptron_per_mouse(mouseNo).percent_correct)
                per_ii=handles_out.perceptron_per_mouse(mouseNo).percent_correct(percent_correct_ii).no_trials;
                
                if per_ii>=20
                    figNo=figNo+1
                    try
                        close(figNo)
                    catch
                    end
                    
                    figure(figNo)
                    
                    hold on
                    
                    perceptron_correct_shuffled=zeros(1,length(t));
                    perceptron_correct_shuffled(1,:)=handles_out.perceptron_per_mouse(mouseNo).percent_correct(percent_correct_ii).perceptron_correct_shuffled;
                    perceptron_correct=zeros(1,length(t));
                    perceptron_correct(1,:)=handles_out.perceptron_per_mouse(mouseNo).percent_correct(percent_correct_ii).perceptron_correct;
                    
                    
                    per95=prctile(perceptron_correct_shuffled,95);
                    per5=prctile(perceptron_correct_shuffled,5);
                    CIsh=[mean(perceptron_correct_shuffled)-per5 per95-mean(perceptron_correct_shuffled)]';
                    [hlCR, hpCR] = boundedline([t(1) t(end)],[mean(perceptron_correct_shuffled) mean(perceptron_correct_shuffled)], CIsh', 'r');
                    
                    plot(t',perceptron_correct,'-k')
                    
                    %Odor on markers
                    plot([0 0],[0 100],'-k')
                    odorhl=plot([0 2.5],[20 20],'-k','LineWidth',5);
                    plot([2.5 2.5],[0 100],'-k')
                    
                    title(['Percent correct prediction by perceptron for mouse No ' num2str(mouseNo) ' ' handles_out.drgbchoices.per_lab{percent_correct_ii}])
                    xlabel('Time (sec)')
                    ylabel('Percent correct')
                else
                    fprintf(1, ['Perceptron not processed for mouse No %d ' handles_out.drgbchoices.per_lab{percent_correct_ii} ' because there were only %d trials (fewer than 20 trials)\n'],mouseNo,per_ii);
                end
            end
        end
    end
    
end
function data_this_trial=drgGetThisTrialEDF(handles,trialNo)
ptr_file=matfile([handles.drg.drta_p.fullName(1:end-4) '_edf.mat']);

data_this_trial=zeros(handles.drg.draq_p.sec_per_trigger*handles.drg.draq_p.ActualRate,22);
data_this_trial(:,1:ptr_file.no_columns)=ptr_file.data_out(handles.drg.draq_d.trial_ii_start(trialNo):handles.drg.draq_d.trial_ii_end(trialNo),:);





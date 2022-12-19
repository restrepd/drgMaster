function data_this_trial=drgGetThisTrialEDF(handles,trialNo)

if handles.upload_edf==1
    data_this_trial=zeros(handles.drg.draq_p.sec_per_trigger*handles.drg.draq_p.ActualRate,22);
    data_this_trial(:,1:4)=handles.data_out(handles.drg.draq_d.trial_ii_start(trialNo):handles.drg.draq_d.trial_ii_end(trialNo),:);
else
    data_this_trial=zeros(handles.drg.draq_p.sec_per_trigger*handles.drg.draq_p.ActualRate,22);
    data_this_trial(:,1:handles.ptr_file.no_columns)=handles.ptr_file.data_out(handles.drg.draq_d.trial_ii_start(trialNo):handles.drg.draq_d.trial_ii_end(trialNo),:);
    pfft=1;
end




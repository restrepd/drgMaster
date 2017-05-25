function data_this_trial = drgGetThisTrial(handles,evNo)
%This function gets the data for all channels for this trial

sessionNo=handles.sessionNo;

trialNo=find(handles.drg.session(sessionNo).trial_start<handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),1,'last');

switch handles.drg.session(sessionNo).draq_p.dgordra
    
    case 1
        %This is dra, this has not been updated, and will not work
        offset=handles.draq_p.no_spike_ch*2*(sum(handles.draq_d.samplesPerTrial(1:handles.p.trialNo))-handles.draq_d.samplesPerTrial(handles.p.trialNo));
        fseek(fid,offset,'bof');
        data_vec=fread(fid,handles.draq_p.no_spike_ch*handles.draq_d.samplesPerTrial(handles.p.trialNo),'uint16');
        szdv=size(data_vec);
        data_this_trial=reshape(data_vec,szdv(1)/handles.draq_p.no_spike_ch,handles.draq_p.no_spike_ch);
    case 2
        %This is dg
        if handles.read_entire_file==1
            data_this_trial_vec=handles.data_dg(floor(handles.drg.session(sessionNo).draq_p.ActualRate*handles.drg.session(sessionNo).draq_p.sec_per_trigger...
                *handles.drg.session(sessionNo).draq_p.no_chans*(trialNo-1)+1):...
                floor(handles.drg.session(sessionNo).draq_p.ActualRate*handles.drg.session(sessionNo).draq_p.sec_per_trigger*...
                handles.drg.session(sessionNo).draq_p.no_chans*trialNo)-2000);
        else
            fid=fopen(handles.drg.drta_p.fullName,'rb');
            bytes_per_native=2;     %Note: Native is unit16
            %Note: Since 2013 DT3010 is acquiring at a rate that is not an integer.
            %However, the program saves using an integer number of bytes! VERY
            %important to use uint64(handles.draq_p.ActualRate) here!!! Nick George
            %solved this problem!
            size_per_ch_bytes=handles.drg.session(sessionNo).draq_p.sec_per_trigger*uint64(handles.drg.session(sessionNo).draq_p.ActualRate*bytes_per_native);
            no_unit16_per_ch=size_per_ch_bytes/bytes_per_native;
            trial_offset=handles.drg.session(sessionNo).draq_p.no_chans*size_per_ch_bytes*(trialNo-1);
            status=fseek(fid, trial_offset, 'bof');
            data_this_trial_vec=fread(fid,no_unit16_per_ch*handles.drg.session(sessionNo).draq_p.no_chans,'uint16');
            fclose(fid); 
        end
        
        data_this_trial=zeros(floor(handles.drg.session(sessionNo).draq_p.ActualRate*handles.drg.session(sessionNo).draq_p.sec_per_trigger),handles.drg.session(sessionNo).draq_p.no_chans);
        for ii=1:handles.drg.session(sessionNo).draq_p.no_chans
            data_this_trial(1:end-2000,ii)=data_this_trial_vec(floor((ii-1)*handles.drg.session(sessionNo).draq_p.ActualRate*handles.drg.session(sessionNo).draq_p.sec_per_trigger)+1:...
                floor((ii-1)*handles.drg.session(sessionNo).draq_p.ActualRate*handles.drg.session(sessionNo).draq_p.sec_per_trigger)...
                +floor(handles.drg.session(sessionNo).draq_p.ActualRate*handles.drg.session(sessionNo).draq_p.sec_per_trigger)-2000);
        end
        
    case 3
        data_this_trial=drgGetThisTrialRHD(handles,trialNo);
end






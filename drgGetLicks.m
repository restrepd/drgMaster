function [these_events]=drgGetLicks(Fs,threshold,lick_voltage)
%Get licks

the_end=0;

these_events=[];
no_evs_this_trial=0;

this_lick_ii=0;
these_lick_times=[];

inter_lick_intervals=[];
iii_ili=0;

smallest_inter_lick_interval=0.02;  %Note: this is used to reject "lick" events due to noise

ii=1;

%Find the events (licks)
while the_end==0
    next_event=find(lick_voltage(ii:end)>threshold,1,'first');
    if isempty(next_event)
        the_end=1;
    else
        
        ii=ii+next_event-1;
        
        %Exclude if the inter event interval is too
        %small due to noise in the lick signal
        
        %Find the inter lick interval
        this_lick_ii=this_lick_ii+1;
        these_lick_times(this_lick_ii)=(ii/Fs);
        if this_lick_ii>1
            %record the inter lick interval
            iii_ili=iii_ili+1;
            inter_lick_intervals(iii_ili)=these_lick_times(this_lick_ii)-these_lick_times(this_lick_ii-1);
        else
            iii_ili=iii_ili+1;
            inter_lick_intervals(iii_ili)=these_lick_times(this_lick_ii);
        end
        
        %Enter the event (lick) in the timecourse only if it is
        %not within a burst of high frequency noise
        
        if inter_lick_intervals(iii_ili)>smallest_inter_lick_interval
            no_evs_this_trial=no_evs_this_trial+1;
            these_events(no_evs_this_trial)=ii/Fs;
        end
        
        
        end_event=find(lick_voltage(ii:end)<threshold,1,'first');
        if isempty(end_event)
            the_end=1;
        else
            ii=ii+end_event-1;
        end
    end
end

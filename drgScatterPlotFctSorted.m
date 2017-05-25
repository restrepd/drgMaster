function drgScatterPlotFctSorted(handles)
%Scatter plot


textout='drgScatterPlotFctSorted'
tic

sessionNo=handles.drg.unit(handles.unitNo).sessionNo;


try
    close 1
catch
end

evTypeNo=handles.evTypeNo;
unitNo=handles.unitNo;
firstTr=handles.trialNo;
lastTr=handles.lastTrialNo;
drg=handles.drg;


hFig=figure(1);
set(hFig, 'units','normalized','position',[.025 .1 .4 .8])

szblocks=size(handles.drg.session(handles.drg.unit(handles.unitNo).sessionNo).blocks);
numEvs=zeros(1,szblocks(1));
blkindx=1:szblocks(1);

spikes1=handles.drg.unit(handles.unitNo).spike_times;
szindx1=size(spikes1);




num_trials=0;
for trNo=firstTr:lastTr
    
    evNo = drgFindEvNo(handles,trNo,sessionNo,evTypeNo);
    if evNo~=-1
        excludeTrial=drgExcludeTrial(drg,drg.unit(unitNo).channel,drg.session(sessionNo).events(evTypeNo).times(evNo),sessionNo);
        
        if excludeTrial==0
            num_trials=num_trials+1;
            BFR(num_trials)=-sum((spikes1>handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo)+handles.time_start+handles.time_pad)&...
                (spikes1<=handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo)))/(handles.time_start+handles.time_pad);
            event_no(num_trials)=evNo;
        end
    end
end

to_sort=[BFR' event_no'];
sorted_BFR=sortrows(to_sort,1);
sorted_evNo=sorted_BFR(:,2);

hold on


plot_y=szblocks(1)*10+num_trials*5+10;

previous_block=-1;

% ylim([0 plot_y*1.01]);
y_min=200000;
y_max=-20000;
xlim([handles.time_start+handles.time_pad handles.time_end-handles.time_pad]);



how_many=0;


num_plotted=0;

for ii=1:length(sorted_evNo)
    
         evNo=sorted_evNo(ii);
            
            perieventindx=[];
            perieventSpikes=[];
            perieventindx=(spikes1>(handles.drg.session(handles.drg.unit(handles.unitNo).sessionNo).events(handles.evTypeNo).times(evNo)+handles.time_start+handles.time_pad))&(spikes1<(handles.drg.session(handles.drg.unit(handles.unitNo).sessionNo).events(handles.evTypeNo).times(evNo)+handles.time_end-handles.time_pad));
            perieventSpikes=spikes1(perieventindx)-handles.drg.session(handles.drg.unit(handles.unitNo).sessionNo).events(handles.evTypeNo).times(evNo);
            szpcs=size(perieventSpikes);
            
            
            if szpcs(2)>0
                for ii=1:szpcs(2)
                    %if ((perieventSpikes(ii)>-x_from)&(perieventSpikes(ii)<x_to))
                    t=[perieventSpikes(ii) perieventSpikes(ii)];
                    y=[plot_y plot_y-5];
                    plot(t,y,'-');
                    if plot_y>y_max
                        y_max=plot_y;
                    end
                    if plot_y-5<y_min
                        y_min=plot_y-5;
                    end
                    %end
                end
                plot_y=plot_y-5;
                
                
                num_plotted=num_plotted+1;
                if num_plotted>length(handles.drg.session(handles.drg.unit(handles.unitNo).sessionNo).events(handles.evTypeNo).times)
                    num_plotted=noTr+1;
                end
            end
       
    
    
    if y_min<y_max
        ylim([y_min-0.05*(y_max-y_min) y_max+0.05*(y_max-y_min)])
    end
    
end %for evNo

title(['Scatterplot for ' handles.drg.session.eventlabels{handles.evTypeNo}])
xlabel('Time (sec)')
text(-2,y_min-0.03*(y_max-y_min),'last')
text(-2,y_max+0.02*(y_max-y_min),'first')

toc
function drgDisplayBatchSync(handles)

%This function displays the LFP power spectrum for drgRunBatch-generated
%data



% 1 Show the difference between events for each group/ percent correct
%  bin in a separate graph
%
% 2 Show the difference between groups for each event/percent correct
%  bin in a separate graph

%Which time window do you want displayed
winNo=1;


close all
warning('off')


%Ask user for the drgb output .mat file and load those data
[handles.drgb.outFileName,handles.PathName] = uigetfile('*.mat','Select the drgb output file');
load([handles.PathName handles.drgb.outFileName])


%Initialize the variables

%Determine how many LFPs are included for each evTypeNo, time window, group, etc
p_values=[];
percent_synch=[];
MI_Tort=[];
groupNo=[];
skew=[];
for unit_pair=1:handles_drgb.drgb.corrpair_no
   if (handles_drgb.drgb.corrpair(unit_pair).ref_u_single==1)&(handles_drgb.drgb.corrpair(unit_pair).partner_u_single==1)
       %Both are single units
       p_values=[p_values handles_drgb.drgb.corrpair(unit_pair).p_value];
       percent_synch=[percent_synch handles_drgb.drgb.corrpair(unit_pair).percent_synch];
       MI_Tort=[MI_Tort handles_drgb.drgb.corrpair(unit_pair).MI_Tort];
       groupNo=[groupNo handles_drgb.drgb.corrpair(unit_pair).groupNo];
       skew=[skew handles_drgb.drgb.corrpair(unit_pair).skewness];
   end    
end
pFDR=drsFDRpval(p_values);

%Perform per unit analysis
basalFR=[];
groupNoU=[];
SingleUnit=[];
SingleUnit_Fee=[];
perViol=[];
for unitNo=1:handles_drgb.drgb.unit_no
       basalFR=[basalFR handles_drgb.drgb.unit(unitNo).basalFR];
       groupNoU=[groupNoU handles_drgb.drgb.unit(unitNo).groupNo];
       SingleUnit=[SingleUnit handles_drgb.drgb.unit(unitNo).SingleUnit];
       SingleUnit_Fee=[SingleUnit_Fee handles_drgb.drgb.unit(unitNo).SingleUnit_Fee];
       perViol=[perViol handles_drgb.drgb.unit(unitNo).perViol];
end

%These are the colors for the different lines
these_colors{1}='b';
these_colors{2}='r';
these_colors{3}='m';
these_colors{4}='g';
these_colors{5}='y';
these_colors{6}='k';
these_colors{7}='c';
these_colors{8}='k';

these_lines{1}='-b';
these_lines{2}='-r';
these_lines{3}='-m';
these_lines{4}='-g';
these_lines{5}='-y';
these_lines{6}='-k';
these_lines{7}='-c';
these_lines{8}='-k';

these_circles{1}='*b';
these_circles{2}='*r';
these_circles{3}='*m';
these_circles{4}='*g';
these_circles{5}='*y';
these_circles{6}='*k';
these_circles{7}='*c';
these_circles{8}='*k';


%basalFR
figNo=1;
try
    close(figNo)
catch
end
figure(figNo)

[f_bFR_WT,x_bFR_WT] = drg_ecdf(basalFR((groupNoU==1)));
plot(x_bFR_WT,f_bFR_WT,these_lines{1})

hold on

[f_bFR_Null,x_bFR_Null] = drg_ecdf(basalFR((groupNoU==2)));
plot(x_bFR_Null,f_bFR_Null,these_lines{2})

%Now do the figures and statistical analysis for synchrony

%MI Tort
figNo=figNo+1;
try
    close(figNo)
catch
end
figure(figNo)

[f_mi_WT,x_mi_WT] = drg_ecdf(MI_Tort((groupNo==1)&(p_values<=pFDR)&(~isnan(MI_Tort))&(~isinf(MI_Tort))));
plot(x_mi_WT,f_mi_WT,these_lines{1})

hold on

[f_mi_Null,x_mi_Null] = drg_ecdf(MI_Tort((groupNo==2)&(p_values<=pFDR)&(~isnan(MI_Tort))&(~isinf(MI_Tort))));
plot(x_mi_Null,f_mi_Null,these_lines{2})


%percent synch
figNo=figNo+1;
try
    close(figNo)
catch
end
figure(figNo)

[f_ps_WT,x_ps_WT] = drg_ecdf(percent_synch((groupNo==1)&(p_values<=pFDR)));
plot(x_ps_WT,f_ps_WT,these_lines{1})

hold on

[f_ps_Null,x_ps_Null] = drg_ecdf(percent_synch((groupNo==2)&(p_values<=pFDR)));
plot(x_ps_Null,f_ps_Null,these_lines{2})

%skewness
figNo=figNo+1;
try
    close(figNo)
catch
end
figure(figNo)

[f_sk_WT,x_sk_WT] = drg_ecdf(skew((groupNo==1)&(p_values<=pFDR)));
plot(x_sk_WT,f_sk_WT,these_lines{1})

hold on

[f_sk_Null,x_sk_Null] = drg_ecdf(skew((groupNo==2)&(p_values<=pFDR)));
plot(x_sk_Null,f_sk_Null,these_lines{2})

%Percent of synchronous unit pairs
perc_synch_up_WT=100*sum((groupNo==1)&(p_values<=pFDR))/sum(groupNo==1)

perc_synch_up_Null=100*sum((groupNo==2)&(p_values<=pFDR))/sum(groupNo==2)

pffft=1


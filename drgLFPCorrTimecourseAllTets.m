function handles=drgLFPCorrTimecourseAllTets(handles)

handles.displayData=0;
handles.corr_out=[];

gcp

%handles.burstLFPNo reference LFP (hippocampus)
%handles.peakLFPNo the other LFP

%Between tetrodes
handles.corr_out.no_between_pairs=0;
par_out_LFPcorr=[];


for ref_tetNo=1:4
    for other_tetNo=ref_tetNo+1:4
        for ref_tet_electNo=4*(ref_tetNo-1)+1:4*(ref_tetNo-1)+4
            handles.burstLFPNo=ref_tet_electNo;
            for other_tet_electNo=4*(other_tetNo-1)+1:4*(other_tetNo-1)+4
                handles.peakLFPNo=other_tet_electNo;
                handles.corr_out.no_between_pairs=handles.corr_out.no_between_pairs+1;
                handles.corr_out.between.burstLFPNo(handles.corr_out.no_between_pairs)=handles.burstLFPNo;
                handles.corr_out.between.peakLFPNo(handles.corr_out.no_between_pairs)=handles.peakLFPNo;
                handles.corr_out.between.refTetNo(handles.corr_out.no_between_pairs)=ref_tetNo;
                handles.corr_out.between.otherTetNo(handles.corr_out.no_between_pairs)=other_tetNo;
                par_out_LFPcorr(handles.corr_out.no_between_pairs).LFPcorr=[];
                
            end
        end
    end
end

%Run drgLFPCorrTimecourse
no_between_pairs=handles.corr_out.no_between_pairs;
parfor ii_between_pairs=1:no_between_pairs
% for ii_between_pairs=1:no_between_pairs
    handlespf=struct();
    handlespf=handles;
    handlespf.burstLFPNo=handlespf.corr_out.between.burstLFPNo(ii_between_pairs);
    handlespf.peakLFPNo=handlespf.corr_out.between.peakLFPNo(ii_between_pairs);
    handlespf=drgLFPCorrTimecourse(handlespf);
    par_out_LFPcorr(ii_between_pairs).LFPcorr=handlespf.drgb.LFPcorr;
    fprintf(1, ['Between tetrodes: Reference electrode %d, other electrode %d\n'], handlespf.burstLFPNo,handlespf.peakLFPNo);
end

%Save the output to handles
for ii_between_pairs=1:no_between_pairs
    handles.corr_out.between.LFPcorr(ii_between_pairs).LFPcorr=par_out_LFPcorr(ii_between_pairs).LFPcorr;
end

save('tempfile.mat','handles','-v7.3')

%Within tetrodes
handles.corr_out.no_within_pairs=0;
par_out_LFPcorr=[];


for ref_tetNo=1:4
    
        for ref_tet_electNo1=4*(ref_tetNo-1)+1:4*(ref_tetNo-1)+4
            handles.burstLFPNo=ref_tet_electNo1;
            for ref_tet_electNo2=ref_tet_electNo1+1:4*(ref_tetNo-1)+4
                handles.peakLFPNo=ref_tet_electNo2;
                handles.corr_out.no_within_pairs=handles.corr_out.no_within_pairs+1;
                handles.corr_out.within.burstLFPNo(handles.corr_out.no_within_pairs)=handles.burstLFPNo;
                handles.corr_out.within.peakLFPNo(handles.corr_out.no_within_pairs)=handles.peakLFPNo;
                handles.corr_out.within.refTetNo(handles.corr_out.no_within_pairs)=ref_tetNo;
                par_out_LFPcorr(handles.corr_out.no_within_pairs).LFPcorr=[];
                
            end
        end
    
end

%Run drgLFPCorrTimecourse
no_within_pairs=handles.corr_out.no_within_pairs;
parfor ii_within_pairs=1:no_within_pairs
% for ii_within_pairs=1:no_within_pairs
    handlespf=struct();
    handlespf=handles;
    handlespf.burstLFPNo=handlespf.corr_out.within.burstLFPNo(ii_within_pairs);
    handlespf.peakLFPNo=handlespf.corr_out.within.peakLFPNo(ii_within_pairs);
    handlespf=drgLFPCorrTimecourse(handlespf);
    par_out_LFPcorr(ii_within_pairs).LFPcorr=handlespf.drgb.LFPcorr;
    fprintf(1, ['Within tetrodes: Reference electrode %d, other electrode %d\n'], handlespf.burstLFPNo,handlespf.peakLFPNo);
end

%Save the output to handles
for ii_within_pairs=1:no_within_pairs
    handles.corr_out.within.LFPcorr(ii_within_pairs).LFPcorr=par_out_LFPcorr(ii_within_pairs).LFPcorr;
end

save('tempfile.mat','handles','-v7.3')

%Now plot the max_rho_t_lag
try
    close 3
catch
end


hFig3 = figure(3);
set(hFig3, 'units','normalized','position',[.05 .05 .5 .8])

time=handles.corr_out.within.LFPcorr(ii_within_pairs).LFPcorr.time;
sp_ii=0;
for ref_tetNo=1:4
    for other_tetNo=ref_tetNo+1:4

        %         ii=0;
        sp_ii=sp_ii+1;
        subplot(4,2,sp_ii)
        hold on
        for ii_between_pairs=1:no_between_pairs
            if (handles.corr_out.between.refTetNo(ii_between_pairs)==ref_tetNo)&...
                    (handles.corr_out.between.otherTetNo(ii_between_pairs)==other_tetNo)
                %                 ii=ii+1;
                %                 max_rho_t_lag(ii,:)=handles.corr_out.between.LFPcorr(ii_between_pairs).LFPcorr.max_rho_t_lag;
                %                 plot(time,handles.corr_out.between.LFPcorr(ii_between_pairs).LFPcorr.max_rho_t_lag,'-b','LineWidth',0.5)
                max_rho_t_lag=handles.corr_out.between.LFPcorr(ii_between_pairs).LFPcorr.max_rho_t_lag;
                plot(time(abs(max_rho_t_lag)<0.06),max_rho_t_lag(abs(max_rho_t_lag)<0.06),'.b','LineWidth',0.5)
            end
        end
        
        title(['Tetrodes ' num2str(ref_tetNo) ' and ' num2str(other_tetNo)])
        xlabel('Time (sec)')
        ylabel('Lag (sec)')
        
        %         mean_max_rho_t_lag=mean(max_rho_t_lag,1);
        %         tempCI = bootci(1000, {@mean, max_rho_t_lag})';
        %         CI(:,1)=mean_max_rho_t_lag'-tempCI(:,1);
        %         CI(:,2)=tempCI(:,2)-mean_max_rho_t_lag';
        %
        %        [hlCR, hpCR] = boundedline(time,mean_max_rho_t_lag, CI, 'b');
        
    end
end



pffft=1;
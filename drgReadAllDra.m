function data_dra=drgReadAllDra(drafile,sessionNo,drg)
%VERY IMPORTANT: Read all the data from the dg file!!
%reading with an offset will not work

data_dra=[];
fid=fopen(drafile,'r');
no_unit16_per_ch=floor(drg.session(sessionNo).draq_p.sec_per_trigger*drg.session(sessionNo).draq_p.ActualRate);
data_dra=fread(fid,no_unit16_per_ch*drg.session(sessionNo).draq_p.no_spike_ch*drg.session(sessionNo).draq_d.noTrials,'uint16');
fclose(fid)
 

%drgDisplayEDFinfo reads the edfinfo for the .edf file
[FileName,PathName] = uigetfile({'*.edf'},'Select the .edf file with power analysis');
fprintf(1, ['\ndrgDisplayBatchMultiDayLFP run for ' FileName '\n\n']);

einf=edfinfo(handlespf.drg.drta_p.fullName)

   
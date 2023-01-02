 

%drgDisplayEDFinfo reads the edfinfo for the .edf file
[FileName,PathName] = uigetfile({'*.edf'},'Select the .edf file with power analysis');
fprintf(1, ['\ndrgDisplayEDFinfo run for ' FileName '\n\n']);

einf=edfinfo([PathName FileName])

   
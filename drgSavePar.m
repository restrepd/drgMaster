function drgSavePar(path,file,this_lfp_per_file,filNum)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
        fprintf(1, 'Saved temp file for file number: %d\n',filNum);
        if isempty(dir(path(1:end-1)))
           mkdir(path(1:end-1)) 
        end
        save([path 'temp_' file],'this_lfp_per_file','-v7.3')

end


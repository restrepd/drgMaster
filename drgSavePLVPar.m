function drgSavePLVPar(path,file,this_PLV_per_mouse,mouse_no)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
        fprintf(1, 'Saved temp file for mouse number: %d\n',mouse_no);
        if isempty(dir(path(1:end-1)))
           mkdir(path(1:end-1)) 
        end
        save([path 'temp_' file],'this_PLV_per_mouse','-v7.3')

end


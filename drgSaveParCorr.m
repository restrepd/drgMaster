function drgSaveParCorr(path,file,this_LFPcorr)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%         fprintf(1, ['Saved ' file '\n']);
        if isempty(dir(path(1:end-1)))
           mkdir(path(1:end-1)) 
        end
        save([path file],'this_LFPcorr','-v7.3')

end


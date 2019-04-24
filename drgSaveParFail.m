function drgSaveParFail(path,file,filNum,handlespf)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
        fprintf(1, 'Saved failure file for file number: %d\n',filNum);
        if isempty(dir(path(1:end-1)))
           mkdir(path(1:end-1)) 
        end
        save([path 'failed_' file],'filNum','handlespf','-v7.3')

end


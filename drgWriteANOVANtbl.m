function drgWriteANOVANtbl(anovanTbl,fileID)
%Writes the output for ANOVAN to a file

%Which is the longest text descriptor
max_length=0;
for ii=1:size(anovanTbl,1)
    max_length=max([max_length length(anovanTbl{ii,1})]);
end

space_between=15;

%Write the headings
fprintf(fileID, ['\n']);

fprintf(fileID, anovanTbl{1,1});
for jj=1:((max_length+space_between)-length(anovanTbl{1,1}))
    fprintf(fileID, ' ');
end

for column_ii=2:size(anovanTbl,2)
    fprintf(fileID, anovanTbl{1,column_ii});
    for jj=1:((space_between)-length(anovanTbl{1,column_ii}))
        fprintf(fileID, ' ');
    end
end

fprintf(fileID, '\n');

for row_ii=2:size(anovanTbl,1)
    fprintf(fileID, anovanTbl{row_ii,1});
    for jj=1:((max_length+space_between)-length(anovanTbl{row_ii,1}))-8
        fprintf(fileID, ' ');
    end
    
    for column_ii=2:size(anovanTbl,2)
        for jj=1:((space_between)-length(num2str(anovanTbl{row_ii,column_ii},'%15.5f')))
            fprintf(fileID, ' ');
        end
        fprintf(fileID, num2str(anovanTbl{row_ii,column_ii},'%15.5f'));
    end
    
    fprintf(fileID, '\n');
end



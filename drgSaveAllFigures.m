%This saves all open figures assuming that they are numbered 1 to n
h =  findobj('type','figure');
n = length(h);
PathName='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/laser/LDA/EAPA/';

mkdir(PathName)
for figNo=1:n
       hFig=figure(figNo);
       savefig(hFig,[PathName 'Figure' num2str(figNo)])
end
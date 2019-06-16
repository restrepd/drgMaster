%This saves all open figures assuming that they are numbered 1 to n
h =  findobj('type','figure');
n = length(h);
PathName='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 6 Decision time-LDA all mouse/IAMO JL AM/IAMO_batch_figures/';

mkdir(PathName)
for figNo=1:n
       hFig=figure(figNo);
       savefig(hFig,[PathName 'Figure' num2str(figNo)])
end
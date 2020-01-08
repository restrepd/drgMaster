function meanPA = drgCircMean(all_PA)
%This is a circular mean suitable for the angle in degrees in the -180 to
%180 degree span
meanPA=(180/pi)*circ_axial(circ_mean(all_PA'*pi/180))';
end


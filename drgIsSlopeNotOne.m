function [p_value_slope,p_value_int,p_value_linear]= drgIsSlopeNotOne(x,y)

%Uses fitlm to determine whether the slope is different from one

%Rotate by -45 degrees
x_shift=x-mean(x);
y_shift=y-mean(y);
[x_rotated, y_rotated]= drgRotateMinus45(x_shift,y_shift);
lm = fitlm(x_rotated,y_rotated,'linear');
p_value_slope=lm.Coefficients{2,4};
p_value_int=lm.Coefficients{1,4};

T=anova(lm,'summary');
F=table2array(T(2,4));
p_value_linear=table2array(T(2,5));


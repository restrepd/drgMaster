function [f_out,x_out]=drg_ecdf(y)
%This ecdf takes care of repeated values so that this cumulative
%histo can be used to queary with rand and does not have repeated
%x values so that one can use interp1


[f_out,x] = ecdf(y);
x_int=x(2:end)-x(1:end-1);
mean_x_int=mean(x_int);

for ii=1:length(x)-1
    if x(ii)==x(ii+1)
        if ii==1
            x_out(ii)=x(ii)-mean_x_int/100;
        else
            x_out(ii)=x(ii)-0.01*(x(ii)-x(ii-1));
        end
    else
        x_out(ii)=x(ii);
    end
end
x_out(length(x))=x(end);
x_out=x_out';

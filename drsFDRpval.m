function pFDR = drsFDRpval(p_values)

% Computes the False Discovery Rate critical significance level 
% using the algorithm in Curran-Everett AJP Reg. Int. Comp. Phys. 
% 279:R1-R8, 2000
%
% This FDR procedure is suitable if
% the individual tests are independent or positively dependent
% (e.g., Gaussian variables that are positively correlated or
% independent).
%
% This function was validated by running the data in page R7 of the
% reference above

% Test of drsFDRpval:
%
% p_vals=[0.723 0.369 0.034 0.009 0.008 0.004 0.003 0.003 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001];
% pFDR=drsFDRpval(p_vals)
% 
% Also, look in the web for fdr_bh.m:
%
% [h, crit_p, adj_p]=fdr_bh(pvals,q,method,report);
%
% That function gives the same result for rejecting the null hypotheses (h) as drsFDRpval
% using the default method. 
%
% However, in fdr_bh the critical p value is defined as the last significant p value, different to
% the definition by Curran-Everett. This does not make a difference in rejecting the 
% null hypotheses using p_value<=pFDR
%
% Use of fdr_bh is useful with alternate method 'dep' when the data are negatively correlated  


if (length(p_values)==0)||(length(p_values)==sum(p_values>0.05))
    pFDR=0.05;
else

    %Start with Bonferroni
    pFDR=0.05/length(p_values);

    szp=size(p_values);
    if szp(1)==1
        p_values=p_values';
    end
    szp=size(p_values);
    p_v_sorted=sortrows(p_values,[1]);
    p_indx=[1:szp(1)];
    p_indx=p_indx';
    szpindx=size(p_indx);
    p_den=szpindx(1)*ones(szpindx(1),1);
    di=0.05*p_indx./p_den;

    %Now see if you can increase the pFDR
    for ii=szp(1):-1:1
        if p_v_sorted(ii)<=di(ii);
            pFDR=di(ii);
            break;
        end
    end
end
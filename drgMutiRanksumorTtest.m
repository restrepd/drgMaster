function [output_data] = drgMutiRanksumorTtest(input_data)
% This function performs pairwise t tests or ranksum for a series of data sets
% A t test is used if the data are normal, otherwise a ranksum is used
pvals=[];
output_data.ii_pairs=0;
ii_for_test=length(input_data);
for ii=1:ii_for_test
    for jj=ii+1:ii_for_test
        output_data.ii_pairs=output_data.ii_pairs+1;
        output_data.ii_in_pair(output_data.ii_pairs)=ii;
        output_data.jj_in_pair(output_data.ii_pairs)=jj;
        [output_data.p(output_data.ii_pairs), output_data.r_or_t(output_data.ii_pairs)]=drg_ranksum_or_ttest(input_data(ii).data,input_data(jj).data);
        pvals=[pvals output_data.p(output_data.ii_pairs)];
    end
end

output_data.pFDR = drsFDRpval(pvals);
fprintf(1, ['\n\npFDR = %d \n\n'],output_data.pFDR)

%Now sort the data
these_ii_pairs=[1:output_data.ii_pairs];
to_sort=[pvals' these_ii_pairs'];
sorted_rows=sortrows(to_sort);
output_data.sorted_ii_pairs=sorted_rows(:,2);

is_first=1;
for jj_pair=1:output_data.ii_pairs
    ii=output_data.ii_in_pair(output_data.sorted_ii_pairs(jj_pair));
    jj=output_data.jj_in_pair(output_data.sorted_ii_pairs(jj_pair));
    p=output_data.p(output_data.sorted_ii_pairs(jj_pair));
    if (p>output_data.pFDR)&(is_first==1)
        is_first=0;
        fprintf(1, ['\np values below are > pFDR\n\n'])
    end
    r_or_t=output_data.r_or_t(output_data.sorted_ii_pairs(jj_pair));
    if r_or_t==0
        fprintf(1, ['p value ranksum for ' input_data(ii).description ' vs ' input_data(jj).description ' =  %d\n'],p)
    else
        fprintf(1, ['p value t-test for ' input_data(ii).description ' vs ' input_data(jj).description ' =  %d\n'],p)
    end
end

fprintf(1, ['\n\n'])
end


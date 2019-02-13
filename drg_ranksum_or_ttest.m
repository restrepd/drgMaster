function [p_val,r_or_t] = drg_ranksum_or_ttest(x,y)
%Do a t test or a ranksum depending on whether the distributions are normal
%r_or_t=0 - ranksum
%r_or_t=1 - pairwise t test
%r_or_t=2 - t test
%   Detailed explanation goes here
      if (length(x)<4)||(length(y)<4)
                %adtest does not work with n<4. In that case go the
                %safe way ranksum
                p_val=ranksum(x,y);
                r_or_t=0;
            else
                if (adtest(x)==1)||(adtest(x)==1)
                    p_val=ranksum(x,y);
                    r_or_t=0;
                else
                    if length(x)~=length(y)
                        [h, p_val]=ttest2(x,y);
                        r_or_t=1;
                    else
                        [h, p_val]=ttest(x,y);
                        r_or_t=2;
                    end
                end
            end
end


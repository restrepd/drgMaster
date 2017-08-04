function [permvec]=permVec(vec, replications)
% Create a new vector that is a permutation of the input vector. The length
% of the output can be multiplied by replications, but by default, it is 1.

if nargin < 2
   replications=1;
end

len=length(vec);

permlen=len*replications;
permvec=zeros(1,permlen);

r=1;
while(r<permlen)
    p=randperm(len);
    for i=1:len
        permvec(r)=vec(p(i));
        r=r+1;
    end
end
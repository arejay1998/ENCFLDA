function [ ld ]=get_ld(LD_adjmat)
index=find(LD_adjmat);
[i,j]=ind2sub(size(LD_adjmat),index);
ld = [i,j];
end
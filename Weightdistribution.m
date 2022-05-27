function [ result ] = Weightdistribution(new_lm_association ,new_md_association)

 A = new_lm_association';
 B = new_md_association;
[nlA,ndA] = size(A);

Rscore1_lnc_dis=zeros(nlA,ndA);

for i=1:nlA
        q=bsxfun(@rdivide,repmat(A(i,:),nlA,1).*A,sum(A));
        W_A(i,1:nlA)=1./sum(A,2).*sum(q,2);
end

result= W_A*B;
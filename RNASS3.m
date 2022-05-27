function [ result ] = RNASS3(new_ld_association )
% ndA:the number of diseases in A
% nlA:the number of lncRNAs in A
 A = new_ld_association;
[nlA,ndA] = size(A);

% Initialize matrix Rsocre1_lnc_dis
Rscore1_lnc_dis=zeros(nlA,ndA);

%calculate the corresponding weight matrix W in A
for i=1:nlA
        q=bsxfun(@rdivide,repmat(A(i,:),nlA,1).*A,sum(A));
        W_A(i,1:nlA)=1./sum(A,2).*sum(q,2);
end
%calculate the level of consistency between the contribution of resource moved in both directions 
%obtain the first level of resource score about lncRNA_disease_associations
result= W_A*A;
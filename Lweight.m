function [ ldwight ] = Lweight( new_lm_association, new_md_association )
%mdw = Weightdistribution(new_lm_association, new_md_association);
lmw = Weightdistribution(new_md_association', new_lm_association')';

[row1, col1]=size(lmw);
%[row2, col2]=size(mdw);
for i = 1:row1
    for j = 1:row1
ldwight(i,j) =dot(lmw(i,:),lmw(j,:))/(norm(lmw(i,:))*norm(lmw(j,:)));
    end
end
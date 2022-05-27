function [ ldwight ] = LDweight( new_lm_association, new_md_association )
mdw = Weightdistribution(new_lm_association, new_md_association);
lmw = Weightdistribution(new_md_association', new_lm_association')';

[row1, col1]=size(lmw);
[row2, col2]=size(mdw);
ldwight = zeros(row1,col2);
for i = 1:row1
    for j = 1:col2
ldwight(i,j) =dot(lmw(i,:),mdw(:,j))/(norm(lmw(i,:))*norm(mdw(:,j)));
    end
end
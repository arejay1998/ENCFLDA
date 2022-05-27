function [ ldwight ] = ldwight( new_lm_association, new_md_association )
lmw = RNASS3(new_lm_association);
mdw = RNASS3(new_md_association);
[row1, col1]=size(lmw);
[row2, col2]=size(mdw);
ldwight = zeros(row1,col2);
for i = 1:row1
    for j = 1:col2
ldwight(i,j) =dot(lmw(i,:),mdw(:,j))/(norm(lmw(i,:))*norm(mdw(:,j)));
    end
end
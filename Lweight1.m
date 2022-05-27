function [ ldwight ] = Lweight1( new_lm_association)
%mdw = Weightdistribution(new_lm_association, new_md_association);
%lmw = Weightdistribution(new_md_association', new_lm_association')';
lm = RNASS3(new_lm_association);
[row1, col1]=size(lm);
%[row2, col2]=size(mdw);
for i = 1:row1
    for j = 1:row1
%ldwight(i,j) =dot(lm(i,:),lm(j,:))/(norm(lm(i,:))*norm(lm(j,:)));
ldwight= pdist(lm,'euclidean');
    end
end
end
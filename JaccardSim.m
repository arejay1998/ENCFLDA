function [ C ] = JaccardSim( a)
[rows, cols]=size(a);
C=zeros(rows, rows);
for i=1:rows
    for j=1:rows
A=[
a(i,:);
a(j,:)];
D=pdist(A,'jaccard');
coefficient=1-D;
C(i,j)=coefficient;
    end
end

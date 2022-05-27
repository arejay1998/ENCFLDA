function overallauc = positiontooverall(new_ld_association,new_ld)
%%
%----------------------------------------------------------------------------
% Note: The implementation of this function refers to the code in Chen et al.'s 
% work (DOI: https://academic.oup.com/bioinformatics/article/33/5/733/2736222).
%----------------------------------------------------------------------------
%%
% A: adjacency matrix for the lncRNA-disease associations
% B records the location of non-zero elements in LD

load globalposition.mat;
[n,m]=size(new_ld_association);
[pp,qq]=size(new_ld);
for i=1:pp
    if globalposition(i)>m*n-pp+1
       globalposition(i)=m*n-pp+1;
    end
end
for k=1:m*n-pp+1 
    tp=0;
    for t=1:pp
        if globalposition(1,t)<=k
           tp=tp+1;
        end
    end
    
    fp=k*pp-tp;
    pre(1,k)=tp/(tp+fp);
    rec(1,k)=tp/pp;
end 
figure(2)
plot(pre,rec)
hold on
clear area;
area(1,1)=pre(1,1)*rec(1,1)/2;
for k=2:m*n-pp+1
    area(1,k)=[rec(1,k-1)+rec(1,k)]*[pre(1,k)-pre(1,k-1)]/2;
end
overallauc1=sum(area);
save overallauc1 overallauc1;
end
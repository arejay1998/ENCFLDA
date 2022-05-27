function [weight] = Weight(X)
[m,n]=size(X);  
a=0;  
for j=1:n  
    for i=1:m  
        a=a+X(i,j)^2;  
    end  
    A(1,j)=sqrt(a);  
    a=0;  
end  
A=repmat(A,m,1);  
weight=X./A;  
function [MD_mat_new] = WKNKN6( MD_mat, MM_mat, DD_mat,a)

[rows,cols]=size(MD_mat);
y_m=zeros(rows,cols);  
y_d=zeros(rows,cols);  

knn_network_m = KNN( MM_mat, rows );  %for miRNA
for i = 1 : rows   
        [sort_m,idx_m]=sort(knn_network_m(i,:),2,'descend');   
        for j = 1 : rows 
            y_m(i,:) =  y_m(i,:)+sort_m(1,j)* MD_mat(idx_m(1,j),:); 
        end                      
          y_m(i,:)=y_m(i,:)/ norm(MD_mat(i,:));
end

knn_network_d = KNN( DD_mat , cols );  %for disease
for i = 1 : cols   
        [sort_d,idx_d]=sort(knn_network_d(i,:),2,'descend');
        for j = 1 : cols
            y_d(:,i) =  y_d(:,i)+sort_d(1,j)* MD_mat(:,idx_d(1,j)); 
        end                      
           % y_d(:,i)=y_d(:,i)/sum_d;   
           y_d(:,i)=y_d(:,i)/ norm(MD_mat(:,i));
end


 for i = 1 : rows
     for j = 1 : cols
         MD_mat_new(i,j)=(a*y_m(i,j) +(1-a)* y_d(i,j))/(norm(MM_mat(i,:))+norm(DD_mat(:,j)));
     end    
 end

end

function [ knn_network ] = KNN( network , k )
    [rows, cols] = size( network );
    network= network-diag(diag(network)); 
    knn_network = zeros(rows, cols);
    [sort_network,idx]=sort(network,2,'descend');
    for i = 1 : rows
        knn_network(i,idx(i,1:k))=sort_network(i,1:k);
    end
end



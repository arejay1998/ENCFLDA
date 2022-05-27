clear;

%% Input useful parameters
% 电影特征维数
m= 4;
% 用户数
nu= 8;
% 电影数
nm= 7;

Y= [4,4,-1,-1,1,1,5,-1;5,5,-1,1,-1,-1,-1,-1;-1,4,1,-1,-1,1,5,4;5,-1,2,5,-1,1,2,-1;1,-1,5,4,5,-1,-1,1;1,-1,5,-1,-1,4,-1,-1;-1,1,-1,5,-1,5,1,-1];
R= ones(nm, nu);
index= find(Y== -1);
R(index)= 0;

% 正则化参数
lambda= 0.1;
% 学习率
alpha= 0.01;

%% Initialize
% Random 
X = randn(nm, m);
Theta = randn(nu, m);

% Not random which can't get right result
% c= 3;
% X = c* ones(nm, m);
% Theta = c* ones(nu, m);

%% Gradient Descent
J = 0;
X_grad = zeros(size(X));
Theta_grad = zeros(size(Theta));

num_iter= 50000;
for k= 1:num_iter

    predict=(X*Theta') .* R;
    J=1/2*sum(sum((predict-Y) .^ 2))+lambda/2*sum(sum(Theta .^ 2))+lambda/2*sum(sum(X .^ 2));
    fprintf('J= %d\n', J);

    for i= 1:nm
        idx= find(R(i,:)==1);
        thetatemp= Theta(idx,:);
        Ytemp= Y(i, idx);
        X_grad(i,:)= (X(i,:)*thetatemp'- Ytemp)* thetatemp+ lambda* X(i, :);
        %X(i,:)= X(i,:)- alpha* X_grad(i,:);
    end
    X= X-alpha* X_grad;

    for i= 1:nu
        idx= find(R(:,i)== 1);
        Xtemp= X(idx,:);
        Ytemp= Y(idx,i);
        Theta_grad(i,:)=(Xtemp*Theta(i,:)'-Ytemp)'*Xtemp+lambda*Theta(i,:);
        %Theta(i,:)= Theta(i,:)- alpha* Theta_grad(i,:);
    end
    Theta= Theta-alpha*Theta_grad;

end

%% predict
Y_pred= X* Theta';

%% calculate the square error
idx= find(R== 1);
s_error_mat= (Y_pred- Y).^2;
s_error= sum(s_error_mat(idx));

%% Code for calculate similarity
% msimi= zeros(7, 1);
% for i= 1:7
%     msimi(i)= norm(X(i,:)- X(1,:));
% end

% msimi= zeros(7, 1);
% for i= 1:7
%     msimi(i)= norm(X(i,:)- X(5,:));
% end

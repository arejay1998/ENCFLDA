function [auc] = roc_3(deci,label_y)
[threshold,ind] = sort(deci,'descend');
roc_y = label_y(ind);
stack_x = cumsum(roc_y == 0)/sum(roc_y == 0);
stack_y = cumsum(roc_y == 1)/sum(roc_y == 1);
auc=sum((stack_x(2:length(roc_y))-stack_x(1:length(roc_y)-1)).*stack_y(2:length(roc_y)));

hold on
plot(stack_x,stack_y);
auc
xlabel('False Positive Rate');
ylabel('True Positive Rate');
end

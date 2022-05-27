function [auc] = LOOCV(MiFunSim,DiPheSim,mi2diNetwork)
%Leave-one-out cross validation is implied on TPGLDA
%lncSim represents lncRNA expression similarity as the Spearman correlation coefficient between the expression profiles of each lncRNA pair.
%disSim represents the semantic similarities among all diseases.
[score_ori] = NCPLDA(lncSim02, disSim02, LD_adjmat);
index = find(1 == LD_adjmat);
for i = 1:length(index)
    i
    A(index(i)) = 0;
    [lncRNA,disease]=ind2sub(size(A),index(i));
    if sum(A(lncRNA,:))==0
        A(lncRNA,:)=Sim_lnc(A,lncSim,lncRNA);
    end
    if sum(A(:,disease))==0
        A(:,disease)=Sim_dis(A,disSim,disease);
    end
    [result]=NCPLDA(lncSim02, disSim02, ld_adjmat_new);
    score_ori(index(i)) = result(index(i));
    A = LD_adjmat;
end
pre_label_score = score_ori(:);
save score_ori;
label_y = LD_adjmat(:);
auc=roc(pre_label_score(:),label_y,'b')
end

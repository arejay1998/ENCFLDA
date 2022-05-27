% Load data set

[union_miRNA,L_in_ml,D_in_md,new_ml,new_md,new_ld,new_ml_association,...
    new_md_association,new_ld_association]=get_new_association;
new_lm_association = new_ml_association';
load('ldwight11');
DATA=[];
for A =0.3
    
%% Retain the original correlation matrix and transpose
interaction = new_ld_association;
save interaction interaction;
%% Calculate the cosine similarity matrix and integrate

[m1 ,disSim] = cosSim( new_md_association );
[lncSim1 ,m2] = cosSim( new_lm_association );
disSim01  = JaccardSim( new_md_association' );
lncSim  = JaccardSim( new_lm_association );
%disSim02  = combineSim2(disSim,disSim01);
%lncSim2  = combineSim2(lncSim1,lncSim);
disSim02  = combineSim2(disSim01,disSim);
lncSim2  = combineSim2(lncSim,lncSim1);
%lncSim2 =Lweight(new_lm_association,new_md_association);
%disSim02 =Dweight(new_lm_association,new_md_association);

%% Calculated scoring matrix
%A=135;
matPredict=WKNKN5( predR1, lncSim1, disSim, A);

%% result
[NCP_rank,NCP_rank_known] =Rank_miRNAs( matPredict, new_ld_association, L_in_ml, D_in_md);


%% Leave a cross-validation
index_1 = find(1 == interaction);

for i = 1:length(index_1)
        i
        load interaction;
        interaction(index_1(i))=0;
    
        [result]=WKNKN5(predR1, lncSim2, disSim02, A);
  
        matPredict(index_1(i)) = result(index_1(i));
               
end
% 
%     pre_label_score = NCP(:);
%     label_y = interaction(:);
%     auc = roc_1(pre_label_score,label_y,'red');
%     auc_all(q) = auc;
%     q = q+1;
% end
% x=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9];
% plot(x,auc_all)
% xlabel('¦Á');
% ylabel('AUC');
% grid on
% set(gca,'YLim',[0.7 1]);
% axis([x auc_all],[0 1 0 1])
    pre_label_score =matPredict(:);
%     save pre_label_score_NCPHLDA_lncSim_disSim pre_label_score;
%     save pre_label_score_onlylncSimanddisSim pre_label_score;
   save pre_label_score_NCPLDA pre_label_score;
    load interaction;
    label_y = interaction(:);
    
    auc = roc_1(pre_label_score,label_y);
    %auc = roc_3(pre_label_score,label_y);
    
    aupr = roc_2(pre_label_score,label_y);
    DATA(end+1)= auc;
    save DATA DATA
end
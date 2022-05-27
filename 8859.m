%% Load data set

[union_miRNA,L_in_ml,D_in_md,new_ml,new_md,new_ld,new_ml_association,...
    new_md_association,new_ld_association]=get_new_association;
new_lm_association = new_ml_association';
load('dsim');
%% Retain the original correlation matrix and transpose
interaction = new_ld_association;
save interaction interaction;
%% Calculate the cosine similarity matrix and integrate
%mSim=miRNASS( new_md_association, dsim ); 
%lncSim = miRNASS( new_ld_association, dsim );
%disSim01  = GSD( new_md_association );
%disSim02  = combineSim2(dsim,disSim01);
%lncSim1 = miRNASS( new_lm_association, mSim );
%gama =0.75;
%lncSim2 =plus(gama*(lncSim),(1-gama)*(lncSim1));  
lncSim2 =Lweight(new_lm_association,new_md_association);
disSim02 =Dweight(new_lm_association,new_md_association);

%% Calculated scoring matrix
%a = 0.4;
A=40;
B=0.1;
K=new_lm_association*new_md_association;
%ld_adjmat_new =LDweight(new_lm_association,new_md_association);
KK=K+new_ld_association;
%ld_adjmat_new=Weight(KK);
%ld_adjmat_new=Weight(K);
ld_adjmat_new = RNASS3(KK);
%matPredict=NCPLDA1(lncSim2, disSim02, ld_adjmat_new,a);
matPredict=WKNKN( ld_adjmat_new, lncSim2, disSim02, A, B);

%matPredict=WKNKN( KK, lncSim2, disSim02, A, B);

%% result
[NCP_rank,NCP_rank_known] =Rank_miRNAs( matPredict, new_ld_association, L_in_ml, D_in_md);


%% Leave a cross-validation
index_1 = find(1 == interaction);

for i = 1:length(index_1)
        i
        load interaction;
        interaction(index_1(i))=0;
        %mSim=miRNASS( new_md_association, dsim ); 
        %lncSim = miRNASS( interaction, dsim );
        %disSim01  = GSD( new_md_association );
        %disSim02  = combineSim2(dsim,disSim01);
        %lncSim1 = miRNASS( new_lm_association, mSim );
        %gama =0.75;
        %lncSim2 =plus(gama*(lncSim),(1-gama)*(lncSim1)); 
        %K=new_lm_association*new_md_association;
         KK=K+interaction;
         ld_adjmat_new = RNASS3(KK);
        %ld_adjmat_new=Weight(KK);
        %[result]=NCPLDA1(lncSim2, disSim02, ld_adjmat_new,a);
        [result]=WKNKN( ld_adjmat_new, lncSim2, disSim02, A, B);
        %[result]=WKNKN( KK, lncSim2, disSim02, A, B);

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
    %aupr = roc_2(pre_label_score,label_y);

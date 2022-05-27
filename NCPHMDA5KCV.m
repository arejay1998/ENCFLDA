%% Load data set
[union_miRNA,L_in_ml,D_in_md,new_ml,new_md,new_ld,new_ml_association,...
    new_md_association,new_ld_association]=get_new_association;
new_lm_association = new_ml_association';
load('arejay789');
%% Retain the original correlation matrix and transpose
interaction = new_ld_association;
save interaction interaction;
%% Calculate the cosine similarity matrix and integrate
mSim=miRNASS( new_md_association, dsim ); 
lncSim = miRNASS( new_ld_association, dsim );
disSim01  = GSD( new_md_association );
disSim02  = combineSim2(dsim,disSim01);
lncSim1 = miRNASS( new_lm_association, mSim );
gama =0.75;
lncSim2 =plus(gama*(lncSim),(1-gama)*(lncSim1));    

  
%% Calculated scoring matrix
K=new_lm_association*new_md_association;
KK=K+new_ld_association;
ld_adjmat_new=Weight(KK);
matPredict=NCPLDA(lncSim2, disSim02, ld_adjmat_new);

%% result
%allresult(disease,lncRNA,interaction,NCP);

%% Five fold cross validation
index_1 = find(1 == interaction);
auc=zeros(1,100);
pp = length(index_1);
   for i = 1 : 100
    i
    indices = crossvalind('Kfold', pp, 5); %Randomly divide the data sample into 5 parts
    for j = 1:5  %Cycle 5 times, take the i-th part as the test sample and the other two parts as the training samples
       
        index_2 = find(j == indices);
        load interaction;
        interaction(index_1(index_2)) = 0;
        %mSim=miRNASS( new_md_association, dsim ); 
        lncSim = miRNASS( interaction, dsim );
        %disSim01  = GSD( new_md_association );
        %disSim02  = combineSim2(dsim,disSim01);
        %lncSim1 = miRNASS( new_lm_association, mSim );
        %gama =0.75;
        lncSim2 =plus(gama*(lncSim),(1-gama)*(lncSim1)); 
        %K=new_lm_association*new_md_association;
        KK=K+interaction;
        ld_adjmat_new=Weight(KK);
        [result]=NCPLDA(lncSim2, disSim02, ld_adjmat_new);
        matPredict(index_1(index_2)) = result(index_1(index_2)); 
    end
    pre_label_score = matPredict(:);
    save pre_label_score_NCPHLDA_5kcv pre_label_score;
    load interaction;
    label_y = interaction(:);
    %auc(i) = roc_1(pre_label_score,label_y);
    %auc(i) = roc_3(pre_label_score,label_y);
    auc(i) = roc_2(pre_label_score,label_y);
   end
   auc_avg = mean(auc);
   std(auc);
 
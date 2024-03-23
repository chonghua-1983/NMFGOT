% implement lasso and stabalitity selection
load('processed_data/ERAWIJANTARI_GASTRIC_CANCER_2020/data1.mat')
proData1 = filter_feas(X1,smps1,spe_names, 3, 0, 0); 
proData2 = filter_feas(X2,smps1,mtb_names, 2, 0, 0); 
X1 = proData1.data; microbes = proData1.Features;
X2 = proData2.data; metabolites = proData2.Features;
X2 = log(X2+1); X2 = X2./repmat(sum(X2,2),1,size(X2,2));

% response variable: metabolites; predictor: microbial taxa
load('results/enrichment_analysis/GC/NMFGOT_W1&W2.mat')
[~,I] = sort(W1,1,'descend'); taxa = X1(:,I(1:100,2)); 
[~,I] = sort(W2,1,'descend'); mets = X2(:,I(1:100,2)); % 50 metabolites
taxa_NMFGOT = importdata('results/enrichment_analysis/GC/factor_b2.csv');
taxa_NMFGOT = string(taxa_NMFGOT); taxa_NMFGOT(1) = [];
mtb_NMFGOT = importdata('results/enrichment_analysis/GC/mtb_factor_b2.csv');
mtb_NMFGOT = string(mtb_NMFGOT); mtb_NMFGOT(1) = [];
predictors = array2table(taxa, 'RowNames', smps1, 'VariableNames', taxa_NMFGOT);
response = array2table(mets, 'RowNames', smps1, 'VariableNames', mtb_NMFGOT);
writetable(predictors, 'results/enrichment_analysis/GC/microbes_Lasso.csv','WriteRowNames',true);
writetable(response, 'results/enrichment_analysis/GC/metabolites_Lasso.csv','WriteRowNames',true);

% compute spearman correlation
[corr_spearman, pval_s] = corr(taxa, mets, 'type', 'Spearman');
[corr_pearson, pval_p] = corr(taxa, mets, 'type', 'Pearson');
coef_spearman = corr_spearman .*(pval_s < 0.05);
coef_pearson = corr_pearson .* (pval_p < 0.05);

% import overlap_lasso_stgabsel_ERSD.csv
taxa_tmp = strrep(taxa_NMFGOT, '-', '.');
[~,ind1,ind2] = intersect(overlaplassostabselGC(:,2),taxa_tmp, 'stable');
% manul operation to keep consistence of metabolite names between
% overlaplassostabselGC(:,1) and mtb_NMFGOT
tmp = overlaplassostabselGC(:,1);
[c,id1,id2] = intersect(tmp, mtb_NMFGOT, 'stable');
% prepare for visulzation
stab_spearman = coef_spearman(ind2,id2); stab_pearson = coef_pearson(ind2,id2);
stab_spearman_GC = array2table(stab_spearman, 'RowNames', taxa_tmp(ind2), 'VariableNames', mtb_NMFGOT(id2));
writetable(stab_spearman_GC, 'results/visualization/heatmap/stab_spearman_GC.csv','WriteRowNames',true);

% demo on ERAWIJANTARI_GASTRIC_CANCER_2020 data
% start NMF_opt
load('processed_data/ERAWIJANTARI_GASTRIC_CANCER_2020/data1.mat')
proData1 = filter_feas(X1,smps1,spe_names, 3, 0, 0); 
proData2 = filter_feas(X2,smps1,mtb_names, 2, 0, 0); 
X1 = proData1.data; X2 = proData2.data; 
X2 = log(X2+1);

num_clu = length(unique(label));
num_factor = num_clu;

[W1,H1] = nndsvd(X1',num_factor,0);
[W2,H2] = nndsvd(X2',num_factor,0);
Inits.H1 = H1; Inits.H2 = H2; 
Inits.W1 = W1; Inits.W2 = W2; 

Dist1 = dist2(X1,X1); S1 = affinityMatrix(Dist1,20); 
Dist2 = dist2(X2,X2); S2 = affinityMatrix(Dist2,20); 
S = SNF({S1,S2},20); Inits.S = S; 
clear S W1 H1 W2 H2 

% compute WMD distance between samples
[WMD_D1, ~] = computing_WMD(X1');
[WMD_D2, ~] = computing_WMD(X2');
A1 = affinityMatrix(WMD_D1,20); A1 = A1./max(A1(:));
A2 = affinityMatrix(WMD_D2,20); A2 = A2./max(A2(:));
params = parameter_selection_mfop(X1', X2', A1, A2,num_factor);

[S,W1,W2, H1,H2,obj] = mf_op(X1',X2',A1,A2,Inits,params);
A = Wtrim(S,20); 
[clust,~,~] = getNCluster(A,num_clu,0,3,20); 
if length(unique(clust))== num_clu
    [ac,ari, nmi_value, ~] = CalcMetrics(label, clust); 
else
    [~, clust, ~] = SpectralClustering(A, num_clu); 
    [ac,ari, nmi_value, ~] = CalcMetrics(label, clust); 
end
sil = silhouette_similarityMatrix(A,clust,num_clu);




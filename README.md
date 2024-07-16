# NMFGOT
NMFGOT package for manuscript titled "NMFGOT: A Multi-view Learning Framework for the Microbiome and Metabolome Integrative Analysis with Optimal Transport Plan"

This is an implement of NMFGOT on multi-omics microbiome data.
Running environment：MATLAB R2019b or later.
The external functions used in the manuscript can be found in the folder external/.

1. NMFGOT includes the main functions below: <br>

   "mf_op.m": NMFGOT implementation, integrative analysis for microbiome data and metabolomics data by matrix factorization and optimal transport plan. <br>
   "parameter_selection_mfop.m": parameters selection for NMFGOT. <br>
   "demo1.m": a demo on GASTRIC CANCER(GC) dataset <br>
   "downstream.m": an tutorial for implementing lasso and stabalitity selection to identify the significantly and stable microbe-metabolite associations. <br>
   "tutorial.m": a step-by-step tutorial for implementing NMFGOT, including visualization, parameter selection and so on.<br>
   "external/silhouette_similarityMatrix.m": function for computing silhouette scores metric. <br>
   "external/computing_WMD.m": computing optimal transport distance between samples. <br>
   "data/data1.mat": GC data used in the manuscript to demonstrate the flowchart of NMFGOT. <br>
   "filter_feas.m": filtering features(samples) occurring in samples(features) less than given threshold. <br>
   "HVGs.m": select high variable features based on its average expression and Fano factor
   
2. Reference

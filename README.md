# PreliminaryR21
About: Preliminary data analysis originally for R21, changed to an R01, as a project in collaboration with Dr Connie Cortes lab, looking at molecular exercise signatures influencing prevention and treatment of Alzheimer's disease (AD). 
 
Purpose: Experiments for AD are often performed in mice, which do not perfectly recapitulate disease and response in humans. To improve human-disease interpretation of mouse RNAseq data, transfer learning methods can be utilized to project mouse patterns onto human data.

Data: 
- Mouse dataset consisted of RNAseq from quadricep skeletal muscle with 16 total samples (half on control diet, half on high fat diet, and half the samples were pre exercise and half were post exercise). GEO series GSE97718. 
- Human dataset included RNAseq from vastus lateralis muscle biopsies with a total of 58 samples (30 men, aged 19-30 years, were recruited and divided into two groups: lean (BMI<25, 18.5- 24.1 kg/m2, n=15) and Ov/Ob (BMI≥25, 25.5- 36.9 kg/m2, n=15)). Half the biopsies were collected before exercise, the other half collected post exercise. GEO series GSE108643.

Cortes_ProjectR.Rmd: includes analysis using the ProjectR transfer learning package with the methods: PCA, NMF, and clustering. 
NMF patterns were further used for functional enrichment analysis (FEA) using the gprofiler package.

Cortes_SignatureSearch.Rmd: consists of drug repurposing analysis using gene expression signatures from human pre-post exercise differentially expressed genes. SignatureSearch package was used to explore drug candidates using CMAP, LINCS, gCMAP, and Pearson correlation for querying and reference databases.

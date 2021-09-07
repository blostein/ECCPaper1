#  Evaluating the ecological hypothesis: Early life salivary microbiome assembly predicts dental caries in a longitudinal case-control study
This is the code repository for the paper entitled "Evaluating the ecological hypothesis: Early life salivary microbiome assembly predicts dental caries in a longitudinal case-control study". 

# To run the pipeline in order

--1a) Run all the code in the DADA2 folders to get asv table + taxa tables + representative asv sequences--

1) Run Preprocess/PreProcessMaster.R to preprocess 16S amplicon data from dada2 outputs to phyloseq objects

--2a) Run DMM_Flux.R to get Dirichlet Multinomial model community state type fits across different k values--

2) Run Data_reduction/DataReduceMaster.R to take outputs from steps 1) & 2a) and produce 
     - Community state types at k={4, 5, 6}
     - Random forest fits when using transform = {hellinger, clr}   (COHRA data only) 
     - Weighted cooccurrence networks when using transform={hellinger, clr}
  
--3a) Run the code in the WGS/WGSProcessing/ folder to process whole metagenome data, run humann3 and run metasqueeze (COHRA data only)-- 

3) Run WGS/WGSMaster.R to take outputs from 3a and put them into nice formats for analysis and graphing in R & 
    run DeSEQ on case vs control abundances of KEGG orthologs and taxa (COHRA data only) 

Then, any of the .RMD scripts can be run to produce the figures and tables in the paper
1) Results.RMD will produce results section + main figures and tables 
2) Supplement.RMD will produce supplemental section

# Description of code in each folder 

## Preprocess 
The code in this folder will run the preprocessing on 16S amplicon data outputs from running DADA2 on the COHRA and PRJ data. 
  steps include 
  
      1) CreatePhyloSeq.R - take dada2 outputs and make a phyloseq object, including the reference sequences
      
      2) Decontam.R (COHRA only) - remove contaminant ASVs from phyloseq based on frequency (PicoGreen quantification) & prevalence (negative control) methods, 
      3) FilterSamples.R: 
      
          a) remove samples<1000 reads + samples that are not true samples (controls, empty wells etc)
          b) plot read counts in samples by sample type
          c) plot abundance profiles of mocks and negative controls
          
      4) FilterTaxa.R (abundance + prevalence filter) - remove ASVs present in <5% of the sample which represent <5% of the total abundance of any of the samples that they *are* in 
      5) Metadata cleaning (COHRA data only) - clean the visit metadata for the COHRA2 samples
 Code at the end of PreProcessMaster.R, calculate the number of ASVs and number of samples lost at each of these steps

## Data_reduction
The code in this folder will take phyloseq objects from Preprocess step and run:

     1) CreateNetwork.R - make a weighted correlation network for the clr/hellinger transformed asv matrix.  resulting objects include:
          a) number of modules created
          b) intraStats, network stats (centrality, betweenness) for each taxa and what network each taxa belongs to
          c) relative abundance of each network in each sample in long form 
          d) hub taxa 
          
      2) CreateRF.R (COHRA data only) - run randomforest classifier of case/control status
      
      3) CreateCST.R - run community state type analysis
          a) **important** this only works if the code to create the DMM fits (DMM_Flux.R) has been run on both the COHRA2 and PRJ data first or you have been given the cst fit data
          b) iterates over a choice of k csts to get a matrix of samples assigned to the best fit cst at p>0.8
          
      4) calculate diversity metrics (on non-taxa-filtered phyloseq object)
 
## WGS (COHRA only) 
 The code in this folder will take outputs from metasqueeze and humann3 and put them into nicer formats for R. 
 It will also run DeSeq2 on taxa/KEGG abudances (case vs control) 
 
     1) cleanHumannData.R: take .tsvs from humann3 output and create object containining information on taxa abundance, path abundance and path coverage
     
     2) slimSQMObjects.R: take large SQM object and separate out important tables for futher analysis (taxa, function abundance) 
     
     3) runWGSDeSeq.R: Run Deseq2 to identify taxa/keggs that are differentially abundance between cases and controls 
     
## Analysis 
  The code in this folder is just helper codes that get used multiple times in the analysis or are large and interupt the flow of reading the Results/Supplement.Rmd
  
     1) sampleLoss.R - creates figure and table that tracks count of samples that have/are missing metadata; saliva sample data (16S), plaque/saliva sample (WGS)and performs sanity checks on sample counts
     
     2) createTidyGraph.R - take in correlation graph information and manipulates it into a tidy graph object which is nicer to work with 

## Session info from freida's run 
R version 4.0.2 (2020-06-22)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS  10.16

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
 [1] grid      parallel  stats4    stats     graphics  grDevices utils     datasets  methods  
[10] base     

other attached packages:
 [1] msa_1.22.0                  Biostrings_2.58.0           XVector_0.30.0             
 [4] phangorn_2.7.1.1            ape_5.5                     gtable_0.3.0               
 [7] patchwork_1.1.1             ggvenn_0.1.9                igraph_1.2.6               
[10] ggraph_2.0.5                tidygraph_1.2.0             ggjoy_0.4.1                
[13] ggridges_0.5.3              ggtree_2.4.2                ggalluvial_0.12.3          
[16] aricode_1.0.0               mclust_5.4.7                psych_2.1.6                
[19] compareGroups_4.5.1         openxlsx_4.2.4              readr_1.4.0                
[22] forcats_0.5.1               stringi_1.6.2               broom_0.7.8                
[25] purrr_0.3.4                 splitstackshape_1.4.8       DESeq2_1.30.1              
[28] SummarizedExperiment_1.20.0 Biobase_2.50.0              MatrixGenerics_1.2.1       
[31] matrixStats_0.59.0          GenomicRanges_1.42.0        GenomeInfoDb_1.26.7        
[34] SQMtools_0.6.2              microbiome_1.12.0           vegan_2.5-7                
[37] permute_0.9-5               ggpubr_0.4.0                gridExtra_2.3              
[40] ggrepel_0.9.1               ggthemes_4.2.4              WGCNA_1.70-3               
[43] fastcluster_1.2.3           dynamicTreeCut_1.63-1       tidyr_1.1.3                
[46] data.table_1.14.0           plyr_1.8.6                  scales_1.1.1               
[49] xtable_1.8-4                MLeval_0.3                  cowplot_1.1.1              
[52] e1071_1.7-7                 pROC_1.17.0.1               randomForest_4.6-14        
[55] caret_6.0-88                lattice_0.20-44             ggplot2_3.3.5              
[58] tibble_3.1.2                dplyr_1.0.5                 stringr_1.4.0              
[61] phyloseq_1.34.0             DirichletMultinomial_1.32.0 IRanges_2.24.1             
[64] S4Vectors_0.28.1            BiocGenerics_0.36.1        

loaded via a namespace (and not attached):
  [1] ModelMetrics_1.2.2.2   bit64_4.0.5            knitr_1.33             DelayedArray_0.16.3   
  [5] rpart_4.1-15           RCurl_1.98-1.3         doParallel_1.0.16      generics_0.1.0        
  [9] preprocessCore_1.52.1  RSQLite_2.2.7          mice_3.13.0            proxy_0.4-26          
 [13] chron_2.3-56           bit_4.0.4              webshot_0.5.2          xml2_1.3.2            
 [17] lubridate_1.7.10       assertthat_0.2.1       viridis_0.6.1          gower_0.2.2           
 [21] xfun_0.24              hms_1.1.0              evaluate_0.14          fansi_0.5.0           
 [25] readxl_1.3.1           DBI_1.1.1              geneplotter_1.68.0     tmvnsim_1.0-2         
 [29] Rsolnp_1.16            htmlwidgets_1.5.3      ellipsis_0.3.2         backports_1.2.1       
 [33] annotate_1.68.0        vctrs_0.3.8            here_1.0.1             abind_1.4-5           
 [37] cachem_1.0.5           withr_2.4.2            ggforce_0.3.3          packrat_0.6.0         
 [41] HardyWeinberg_1.7.2    checkmate_2.0.0        treeio_1.14.4          mnormt_2.0.2          
 [45] svglite_2.0.0          cluster_2.1.2          lazyeval_0.2.2         crayon_1.4.1          
 [49] genefilter_1.72.1      recipes_0.1.16         pkgconfig_2.0.3        labeling_0.4.2        
 [53] tweenr_1.0.2           nlme_3.1-152           nnet_7.3-16            rlang_0.4.11          
 [57] lifecycle_1.0.0        extrafontdb_1.0        cellranger_1.1.0       rprojroot_2.0.2       
 [61] polyclip_1.10-0        flextable_0.6.6        aplot_0.0.6            Matrix_1.3-4          
 [65] carData_3.0-4          Rhdf5lib_1.12.1        base64enc_0.1-3        png_0.1-7             
 [69] viridisLite_0.4.0      bitops_1.0-7           rhdf5filters_1.2.1     blob_1.2.1            
 [73] jpeg_0.1-8.1           rstatix_0.7.0          ggsignif_0.6.2         memoise_2.0.0         
 [77] magrittr_2.0.1         zlibbioc_1.36.0        compiler_4.0.2         kableExtra_1.3.4      
 [81] RColorBrewer_1.1-2     cli_3.0.0              ade4_1.7-17            htmlTable_2.2.1       
 [85] Formula_1.2-4          MASS_7.3-54            mgcv_1.8-36            tidyselect_1.1.1      
 [89] yaml_2.2.1             locfit_1.5-9.4         latticeExtra_0.6-29    fastmatch_1.1-0       
 [93] tools_4.0.2            rio_0.5.27             rstudioapi_0.13        uuid_0.1-4            
 [97] foreach_1.5.1          foreign_0.8-81         prodlim_2019.11.13     farver_2.1.0          
[101] Rtsne_0.15             BiocManager_1.30.16    rvcheck_0.1.8          digest_0.6.27         
[105] lava_1.6.9             quadprog_1.5-8         Rcpp_1.0.7             car_3.0-11            
[109] writexl_1.4.0          httr_1.4.2             gdtools_0.2.3          AnnotationDbi_1.52.0  
[113] rsconnect_0.8.18       colorspace_2.0-2       rvest_1.0.0            XML_3.99-0.6          
[117] truncnorm_1.0-8        splines_4.0.2          tidytree_0.3.4         graphlayouts_0.7.1    
[121] multtest_2.46.0        systemfonts_1.0.2      jsonlite_1.7.2         timeDate_3043.102     
[125] ipred_0.9-11           R6_2.5.0               Hmisc_4.5-0            pillar_1.6.1          
[129] htmltools_0.5.1.1      glue_1.4.2             fastmap_1.1.0          BiocParallel_1.24.1   
[133] class_7.3-19           codetools_0.2-18       utf8_1.2.1             curl_4.3.2            
[137] officer_0.3.18         zip_2.2.0              GO.db_3.12.1           Rttf2pt1_1.3.8        
[141] survival_3.2-11        rmarkdown_2.9          biomformat_1.18.0      munsell_0.5.0         
[145] rhdf5_2.34.0           GenomeInfoDbData_1.2.4 iterators_1.0.13       impute_1.64.0         
[149] haven_2.4.1            reshape2_1.4.4         extrafont_0.17        

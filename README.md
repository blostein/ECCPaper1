# ECCPaper1
This is the code repository for the paper entitled " ". 

To run the pipeline in order

-- 1a) Run all the code in the DADA2 folders to get asv table + taxa tables + representative asv sequences --
1) Run Preprocess/PreProcessMaster.R to preprocess 16S amplicon data from dada2 outputs to phyloseq objects

-- 2a) Run DMM_Flux.R to get Dirichlet Multinomial model community state type fits across different k values--
2) Run Data_reduction/DataReduceMaster.R to take outputs from steps 1) & 2a) and produce 
     - Community state types at k={4, 5, 6}
     - Random forest fits when using transform = {hellinger, clr}   (COHRA data only) 
     - Weighted cooccurrence networks when using transform={hellinger, clr}
  
-- 3a) Run the code in the Metasqueeze folder to process whole metagenome data, run humann3 and run metasqueeze (COHRA data only) -- 
3) Run WGS/WGSMaster.R to take outputs from 3a and put them into nice formats for analysis and graphing in R & 
    run DeSEQ on case vs control abundances of KEGG orthologs and taxa (COHRA data only) 

Then, any of the .RMD scripts can be run to produce the figures and tables in the paper
1) Results.RMD will produce results section + main figures and tables (and also supplemental figures)
2) Supplement.RMD will produce supplemental section


- Preprocess 
The code in this folder will run the preprocessing on 16S amplicon data outputs from running DADA2 on the COHRA and PRJ data. 
  steps include 
      1) decontam (COHRA only), 2) sample filtering, 3) taxa filtering (abundance +prevalence filter) 4) Metadata (COHRA data only) 
- Data_reduction
The code in this folder will take phyloseq objects from Preprocess step and run:
      1) community state typing, 2) weighted cooccurance networks 3) random forest (COHRA data only) at various hyperparameters 
 - WGS (COHRA only) 
 The code in this folder will take outputs from metasqueeze and humann3 and put them into nicer formats for R It will also run DeSeq2 on taxa/KEGG abudances (case vs control) 

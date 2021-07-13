# ECCPaper1
This is the code repository for the paper entitled " ". 

# To run the pipeline in order

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
1) Results.RMD will produce results section + main figures and tables 
2) Supplement.RMD will produce supplemental section

# Description of code in each folder 

- Preprocess 
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
      
- Data_reduction
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
      
 - WGS (COHRA only) 
 The code in this folder will take outputs from metasqueeze and humann3 and put them into nicer formats for R. 
 It will also run DeSeq2 on taxa/KEGG abudances (case vs control) 
     1) cleanHumannData.R: take .tsvs from humann3 output and create object containining information on taxa abundance, path abundance and path coverage
     2) slimSQMObjects.R: take large SQM object and separate out important tables for futher analysis (taxa, function abundance) 
     3)runWGSDeSeq.R: Run Deseq2 to identify taxa/keggs that are differentially abundance between cases and controls 
     
 - Analysis 
  The code in this folder is just helper codes that get used multiple times in the analysis or are large and interupt the flow of reading the Results/Supplement.Rmd
     1) sampleLoss.R - creates figure and table that tracks count of samples that have/are missing metadata; saliva sample data (16S), plaque/saliva sample (WGS)and performs sanity checks on sample counts
     2) createTidyGraph.R - take in correlation graph information and manipulates it into a tidy graph object which is nicer to work with 

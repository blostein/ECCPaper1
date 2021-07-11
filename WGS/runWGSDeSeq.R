#a script to run deseq2 on WGS metasqueeze outputs
#load slimmed SQM
load(file.path(data.out, 'MS', 'kegg_and_taxa_slim.Rdata'))

#loop over plaque and saliva 
functions_list=list('plaque_functions'=plaque_kegg_abund, 'salvia_functions'=saliva_kegg_abund, 
                    'plaque_taxa'=plaque_species_abund, 'saliva_taxa'=saliva_species_abund)
KEGGs_list=list('plaque'=plaque_kegg_paths, 'saliva'=saliva_kegg_paths)
results_list=list('plaque_functions'=data.frame(), 'saliva_functions'=data.frame(), 'plaque_taxa'=data.frame(), 'saliva_taxa'=data.frame())
for(i in c(1, 2, 3, 4)){
  #make new dataframe with samples in same order
  metadata=metag_metashort[match(colnames(functions_list[[i]]), metag_metashort$SampleName), ]
  #change rownames to sample names
  rownames(metadata)=metadata$SampleName
  #check order 
  if(all(rownames(metadata)==colnames(functions_list[[i]]))!=T){stop('Rownames in wrong order')}
  #convert to deseq 2 data
  dds=DESeqDataSetFromMatrix(countData=functions_list[[i]], colData=metadata, design= ~ CaseStatus)
  #filter low abundance KEGGS
  keep=rowSums(counts(dds))>=read_limit
  dds=dds[keep,]
  #refactor case status
  dds$CaseStatus=factor(dds$CaseStatus, levels=c('control', 'case'))
  #run deseq 2
  dds2=DESeq(dds)
  #pull results 
  results=results(dds2, name='CaseStatus_case_vs_control')
  if(i %in% c(1, 2)){
    #add kegg path annotation data
    KEGG_paths=as.data.frame(KEGGs_list[[i]])%>%rownames_to_column(var='KEGG')
    results_df=data.frame(results)%>%rownames_to_column(var='KEGG')%>%left_join(KEGG_paths)
    #split kegg annotations into columns
    colnames(results_df)=c(colnames(results_df)[1:7], 'PathAnnotation')
    r2=cSplit(results_df, 'PathAnnotation', '\\|', fixed=F, drop=F)
    r2=cSplit(r2, 'PathAnnotation_01', ';', fixed=F); r2=cSplit(r2, 'PathAnnotation_02', ';', fixed=F)
  }
  if(i %in% c(3, 4)){
    r2=data.frame(results)%>%rownames_to_column(var='Taxa')
  }
  #store in list 
  results_list[[i]]=r2
}

kname=c(saliva_kegg_names, plaque_kegg_names)

r3=bind_rows(results_list[1:2], .id='SampleType')%>%
  mutate(ColorCode=ifelse(padj>0.05, 'padj>0.05', 
                          paste0('padj<0.05: ', 
                                 as.character(PathAnnotation_01_1))), 
         SampleType=gsub('_functions', '', SampleType))%>%
  left_join(data.frame('kname'=kname, 'KEGG'=names(kname))%>%unique())

selected_genera=c('Abiotrophia', 'Actinomyces', 
                  'Candida',  'Neisseria', 'Porphyromonas',
                  'Prevotella', 'Scardovia',  
                  'Streptococcus', 'Veillonella')

t3=bind_rows(results_list[3:4], .id='SampleType')%>%
    mutate(Genus=ifelse(str_detect(Taxa, 'Unclassified'), 
                        gsub('Unclassified ', '', Taxa), 
                        str_extract(Taxa, '\\w*')), 
           ColorCode=ifelse(padj>0.05, 'padj>0.05',
                        ifelse(Genus %in% selected_genera, 
                            paste0('padj<0.05:', as.character(Genus)), 
                            'padj<0.05')), 
           SampleType=gsub('_taxa', '', SampleType))

save(list=c('results_list', 'r3', 't3'), file=file.path(data.out, 'MS', 'deseqresults.Rdata'))

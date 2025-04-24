library(MatrixEQTL)
cancer=c('BLCA','BRCA','CESC','COAD','ESCA',
         'HNSC','KIRC','KIRP','LAML','LGG',
         'LIHC','LUAD','LUSC','PAAD','PCPG',
         'PRAD','SARC','SKCM','STAD','THCA','UCEC')

for(i in 1:length(cancer)){

  apa.pos=read.table(paste0(cancer[i],'_splpos.txt'),header=T,sep='\t')
  cg.pos=read.table(paste0(cancer[i],'_cgpos.txt'),header=T,sep='\t')
  
  file1=c(paste0(cancer[i],'_meth.txt'))
  file2=c(paste0(cancer[i],'_spl.txt'))
  file3=c(paste0(cancer[i],'_clinical.txt'))

  snps = SlicedData$new();
  snps$fileDelimiter = "\t";     
  snps$fileOmitCharacters = "NA"; 
  snps$fileSkipRows = 1;          
  snps$fileSkipColumns = 1;       
  snps$fileSliceSize = 2000;      
  snps$LoadFile(file1)

  gene = SlicedData$new();
  gene$fileDelimiter = "\t";      
  gene$fileOmitCharacters = "NA";
  gene$fileSkipRows = 1;          
  gene$fileSkipColumns = 1;       
  gene$fileSliceSize = 2000;     
  gene$LoadFile(file2)

  cvrt = SlicedData$new();
  cvrt$fileDelimiter = "\t";      
  cvrt$fileOmitCharacters = "NA"; 
  cvrt$fileSkipRows = 1;          
  cvrt$fileSkipColumns = 1;       
  cvrt$fileSliceSize = 2000;      
  cvrt$LoadFile(file3)

  output_cis=paste0(,cancer[i],'CisSig.txt')
  output_trans=paste0(,cancer[i],'TransSig.txt')


  a=as.numeric(nrow(snps))
  b=as.numeric(nrow(gene))

  me=Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output_trans,
    pvOutputThreshold = 0.05/(a*b),
    useModel = modelLINEAR,
    errorCovariance = numeric(),
    verbose = TRUE,
    output_file_name.cis = output_cis,
    pvOutputThreshold.cis = 0.05,
    snpspos = cg.pos,
    genepos = apa.pos,
    cisDist = 1e6,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE
  )

  dfFull = me$param$dfFull
  
  tstat_cis = me$cis$eqtls$statistic;
  r_cis = tstat_cis / sqrt( dfFull + tstat_cis^2 )
  cis=me$cis$eqtls
  cis$r=r_cis
  
  tstat_trans=me$trans$eqtls$statistic
  r_trans=tstat_trans / sqrt( dfFull + tstat_trans^2 )
  trans=me$trans$eqtls
  trans$r=r_trans
  
  cis_sig=cis[cis$FDR<0.05,]

  write.table(cis_sig,output_cis,sep='\t',row.names=F,quote=F)
  write.table(trans,output_trans,sep='\t',row.names=F,quote=F)

}

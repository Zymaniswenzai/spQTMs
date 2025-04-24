#### Immune-related AS Analysis  ####

# Set working directory
setwd('~/08.Immu/03.data')

# Define paths for various input and output files
wd <- c('~/12.spQTM/07.splicing/',       # Splicing and expression data
        '~/12.spQTM/01.spQTM/',           # spQTM results
        '~/12.spQTM/08.Immu/01.TumorPurity/',  # Tumor purity data
        '~/12.spQTM/08.Immu/04.ImmuSplicing/') # Output path

# Define list of cancer types
cancer <- c('BLCA','BRCA','CESC','ESCA','HNSC','KIRC','KIRP','LGG','LIHC','LAML',
            'LUAD','LUSC','PAAD','PCPG','PRAD','SARC','STAD','SKCM','THCA','UCEC')

# Load immune-related gene list
geneList <- read.table('~/12.spQTM/08.Immu/immune_pathway_TREM2GENE.txt', 
                       sep = '\t', header = TRUE)

# Load required libraries
library(msigdbr)
library(fgsea)
library(tibble)
library(dplyr)
library(ggplot2)
library(org.Hs.eg.db)
library(limma)
library(ppcor)

# Loop through each cancer type
for (i in 1:length(cancer)) {
  # Load splicing and expression data
  spl <- read.table(paste0(wd[1], cancer[i], 'spl.txt'), header = TRUE, check.names = FALSE, sep = '\t')
  exp <- read.table(paste0(wd[1], cancer[i], 'exp.txt'), header = TRUE, check.names = FALSE, sep = '\t')

  # Load spQTM results (Cis and Trans)
  tqtm <- read.table(paste0(wd[2], cancer[i], '_TransSig.txt'), header = TRUE, check.names = FALSE, sep = '\t')
  cqtm <- read.table(paste0(wd[2], cancer[i], '_CisSig.txt'), header = TRUE, check.names = FALSE, sep = '\t')

  # Combine unique spQTM genes
  all_genes <- unique(data.frame(id = c(cqtm$gene, tqtm$gene)))

  # Subset splicing data to spQTM genes
  spl_filtered <- merge(all_genes, spl, by.x = 'id', by.y = 'row.names')
  rownames(spl_filtered) <- spl_filtered$id
  spl_filtered$id <- NULL

  # Transpose splicing and expression data
  cpg <- as.data.frame(t(spl_filtered))
  gene <- as.data.frame(t(exp))

  # Load tumor purity data
  tp <- read.table(paste0(wd[3], cancer[i], '_tumorPurity.txt'), sep = '\t', header = TRUE)
  
  # Match tumor purity data with samples
  sample <- data.frame(id = rownames(gene))
  tp <- merge(sample, tp, by = 'id')

  # Initialize results container
  rres <- c()

  # Loop through each splicing event
  for (m in 1:ncol(spl_filtered)) {
    res <- data.frame(splicing = NA, gene = NA, RS = NA)

    # Calculate partial correlation between splicing and each gene, adjusted for tumor purity
    for (j in 1:ncol(gene)) {
      pcor_result <- pcor.test(cpg[, m], gene[, j], tp$TumorPurity)
      res[j, ] <- list(colnames(cpg)[m], colnames(gene)[j],
                       (-log10(pcor_result$p.value)) * sign(pcor_result$estimate))
    }

    # Rank genes based on RS (correlation strength)
    ranked_genes <- res %>%
      arrange(desc(RS)) %>%
      dplyr::select(gene, RS)

    # Map gene symbols to official NCBI symbols
    ranked_genes$NCBI <- alias2SymbolTable(ranked_genes$gene, species = 'Hs')

    # Convert symbols to ENTREZ IDs
    entrez_map <- select(org.Hs.eg.db, keys = ranked_genes$NCBI, 
                         columns = "ENTREZID", keytype = "SYMBOL")
    data2 <- merge(ranked_genes, entrez_map, by.x = 'NCBI', by.y = 'SYMBOL') %>%
      unique() %>%
      dplyr::select(gene = ENTREZID, RS) %>%
      arrange(desc(RS))

    # Format data for GSEA
    ranks <- deframe(data2)

    # Run GSEA with immune gene sets
    fgseaRes <- clusterProfiler::GSEA(geneList = ranks,
                                      TERM2GENE = geneList,
                                      nPerm = 1000,
                                      minGSSize = 5,
                                      maxGSSize = 1000,
                                      pvalueCutoff = 1)

    fgseaResD <- data.frame(fgseaRes)[, 2:10]
    fgseaResD$splicing_RES <- ifelse(fgseaResD$enrichmentScore > 0,
                                     1 - 2 * fgseaResD$pvalue,
                                     2 * fgseaResD$pvalue - 1)
    fgseaResD$splicing_Event <- rep(colnames(cpg)[m], nrow(fgseaResD))

    # Combine results
    rres <- rbind(rres, fgseaResD)

    print(m)  # Progress tracker
  }

  # Save results for the current cancer type
  write.table(rres,
              paste0(wd[4], cancer[i], '_ImmuSpliceScore.txt'),
              sep = '\t', quote = FALSE, row.names = FALSE)
}
